#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ==============================================================================
# OASpepDB: Step 01 - Automated OAS Raw Repertoire Download & Metadata Extractor
# ==============================================================================

import os
import sys
import io
import gzip
import json
import csv
import time
import argparse
import urllib.request
import urllib.error
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock

# Đảm bảo console Windows hiển thị đúng font UTF-8
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# 16 cột metadata trích xuất từ dòng 1 của file OAS
METADATA_FIELDS = [
    "Run", "Link", "Author", "Species", "BSource", "BType", "Longitudinal",
    "Disease", "Subject", "Age", "Vaccine", "Chain", "Unique sequences",
    "Isotype", "Total sequences", "Filename"
]

# Các cột cần giữ lại sau khi lọc Biological QC trong file gz
TRIMMED_COLUMNS = [
    "cdr3_aa", "fr3_tail", "fwr4_aa", "v_call", "d_call", "j_call",
    "v_identity", "sequence_alignment_aa", "Redundancy"
]

csv_lock = Lock()


def format_time(seconds):
    """Đổi số giây thành định dạng HH:MM:SS hoặc MM:SS"""
    if seconds < 0 or seconds > 86400 * 30:
        return "--:--"
    m, s = divmod(int(seconds), 60)
    h, m = divmod(m, 60)
    return f"{h:02d}h {m:02d}m {s:02d}s" if h > 0 else f"{m:02d}m {s:02d}s"


def is_file_complete(filepath):
    """Kiểm tra file đã được lọc hoàn chỉnh chưa (dòng đầu phải bắt đầu bằng cdr3_aa)"""
    if not os.path.exists(filepath) or os.path.getsize(filepath) < 50:
        return False
    try:
        with gzip.open(filepath, "rt", encoding="utf-8", errors="ignore") as f:
            return f.readline().strip().startswith("cdr3_aa,")
    except Exception:
        return False


def download_file(url, temp_path, max_retries=3):
    """Tải file từ URL, có hỗ trợ resume (Range) nếu bị ngắt giữa chừng"""
    headers = {"User-Agent": "Mozilla/5.0"}
    for attempt in range(1, max_retries + 1):
        try:
            existing_size = os.path.getsize(temp_path) if os.path.exists(temp_path) else 0
            req_headers = dict(headers)
            mode = "wb"

            if existing_size > 0:
                req_headers["Range"] = f"bytes={existing_size}-"
                mode = "ab"

            req = urllib.request.Request(url, headers=req_headers)
            with urllib.request.urlopen(req, timeout=60) as resp:
                if existing_size > 0 and resp.status == 200:
                    mode = "wb"

                with open(temp_path, mode) as f_out:
                    while chunk := resp.read(256 * 1024):
                        f_out.write(chunk)

            if os.path.exists(temp_path) and os.path.getsize(temp_path) > 50:
                return True

        except urllib.error.HTTPError as e:
            if e.code == 416:  # File đã tải xong toàn bộ
                return True
            time.sleep(attempt)
        except Exception:
            time.sleep(attempt * 2)

    return False


def trim_file(raw_path, final_path, filename, meta_csv, existing_meta):
    """
    Đọc file raw, lọc Biological QC và ghi stream trực tiếp sang file đích.
    Ghi kèm metadata vào file CSV tổng.
    """
    raw_size = os.path.getsize(raw_path)
    temp_trimmed = final_path + ".tmp.gz"

    try:
        with gzip.open(raw_path, "rt", encoding="utf-8", errors="ignore") as f_in:
            line1 = f_in.readline().strip()
            if not line1:
                return False, raw_size, 0

            # Lấy thông tin metadata ở dòng 1
            meta_dict = json.loads(next(csv.reader([line1]))[0])

            # Đọc danh sách cột ở dòng 2
            line2 = f_in.readline().strip()
            cols = next(csv.reader([line2]))
            col_idx = {c: i for i, c in enumerate(cols)}

            # Kiểm tra các cột bắt buộc
            req_cols = ["productive", "stop_codon", "vj_in_frame", "fwr3_aa", "cdr3_aa", "fwr4_aa"]
            for c in req_cols:
                if c not in col_idx:
                    return False, raw_size, 0

            idx_p = col_idx["productive"]
            idx_sc = col_idx["stop_codon"]
            idx_vj = col_idx["vj_in_frame"]
            idx_f3 = col_idx["fwr3_aa"]
            idx_c3 = col_idx["cdr3_aa"]
            idx_f4 = col_idx["fwr4_aa"]
            min_cols = max(idx_p, idx_sc, idx_vj, idx_f3, idx_c3, idx_f4) + 1

            idx_vc = col_idx.get("v_call", -1)
            idx_dc = col_idx.get("d_call", -1)
            idx_jc = col_idx.get("j_call", -1)
            idx_vi = col_idx.get("v_identity", -1)
            idx_seq = col_idx.get("sequence_alignment_aa", -1)
            idx_red = col_idx.get("Redundancy", -1)

            # Đảm bảo trường Chain trong metadata luôn có giá trị
            if not meta_dict.get("Chain"):
                meta_dict["Chain"] = "Heavy" if "_Heavy_" in filename else "Light"

            # Mở file ghi trực tiếp dạng stream để không tốn RAM
            clean_name = filename[:-3] if filename.endswith(".gz") else filename
            with open(temp_trimmed, "wb") as raw_out:
                with gzip.GzipFile(filename=clean_name, mode="wb", fileobj=raw_out) as gz_out:
                    tw = io.TextIOWrapper(gz_out, encoding="utf-8", newline="", write_through=False)
                    writer = csv.writer(tw)
                    writer.writerow(TRIMMED_COLUMNS)

                    for row in csv.reader(f_in):
                        if len(row) < min_cols:
                            continue

                        # Điều kiện Biological QC
                        if row[idx_p] != "T" or row[idx_sc] != "F" or row[idx_vj] != "T":
                            continue

                        fwr3, cdr3, fwr4 = row[idx_f3], row[idx_c3], row[idx_f4]
                        len_c3 = len(cdr3)
                        if len(fwr3) < 20 or len(fwr4) < 5 or not (5 <= len_c3 <= 40):
                            continue
                        if any(x in cdr3 for x in "X*-"):
                            continue

                        fr3_tail = fwr3[-25:] if len(fwr3) >= 25 else fwr3
                        v_call = row[idx_vc] if idx_vc != -1 and idx_vc < len(row) else ""
                        d_call = row[idx_dc] if idx_dc != -1 and idx_dc < len(row) else ""
                        j_call = row[idx_jc] if idx_jc != -1 and idx_jc < len(row) else ""

                        v_id_val = "0.0"
                        if idx_vi != -1 and idx_vi < len(row) and row[idx_vi]:
                            try:
                                v_id_val = f"{float(row[idx_vi]):.3f}"
                            except Exception:
                                v_id_val = "0.0"

                        seq_aa = row[idx_seq] if idx_seq != -1 and idx_seq < len(row) else ""
                        redundancy = row[idx_red] if idx_red != -1 and idx_red < len(row) else "1"

                        writer.writerow([
                            cdr3, fr3_tail, fwr4, v_call, d_call, j_call,
                            v_id_val, seq_aa, redundancy
                        ])

                    tw.flush()

        trimmed_size = os.path.getsize(temp_trimmed)
        os.replace(temp_trimmed, final_path)

        # Ghi metadata vào file CSV tổng
        with csv_lock:
            if filename not in existing_meta:
                meta_row = [str(meta_dict.get(k, "")) for k in METADATA_FIELDS[:-1]] + [filename]
                with open(meta_csv, "a", encoding="utf-8", newline="") as f_meta:
                    csv.writer(f_meta).writerow(meta_row)
                existing_meta.add(filename)

        return True, raw_size, trimmed_size

    except Exception:
        if os.path.exists(temp_trimmed):
            try:
                os.remove(temp_trimmed)
            except Exception:
                pass
        return False, raw_size, 0


def process_url(entry, data_dir, meta_csv, existing_meta):
    """Quy trình xử lý một URL: Tải về -> Lọc -> Xóa file raw"""
    filename, url = entry
    final_path = os.path.join(data_dir, filename)

    if is_file_complete(final_path):
        return filename, "SKIP", 0, os.path.getsize(final_path)

    temp_raw = final_path + ".downloading"

    # 1. Tải file raw
    if not download_file(url, temp_raw):
        if os.path.exists(temp_raw):
            try:
                os.remove(temp_raw)
            except Exception:
                pass
        return filename, "FAIL_DOWNLOAD", 0, 0

    # 2. Lọc dữ liệu và nén thẳng ra file kết quả
    ok, raw_sz, trimmed_sz = trim_file(temp_raw, final_path, filename, meta_csv, existing_meta)

    # 3. Xóa ngay file raw để tiết kiệm ổ đĩa
    if os.path.exists(temp_raw):
        try:
            os.remove(temp_raw)
        except Exception:
            pass

    return filename, ("OK" if ok else "FAIL_TRIM"), raw_sz, trimmed_sz


def get_urls_from_script(sh_path):
    """Đọc danh sách URL cần tải từ file script shell"""
    urls = []
    if not os.path.exists(sh_path):
        return urls
    with open(sh_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.strip().split()
            for p in parts:
                if p.startswith("http://") or p.startswith("https://"):
                    filename = p.split("/")[-1].strip()
                    urls.append((filename, p))
                    break
    return urls


def load_metadata_set(meta_csv, data_dir):
    """Đọc danh sách file đã có trong file metadata để tránh ghi trùng lặp"""
    done_set = set()
    if not os.path.exists(meta_csv):
        os.makedirs(os.path.dirname(os.path.abspath(meta_csv)), exist_ok=True)
        with open(meta_csv, "w", encoding="utf-8", newline="") as f:
            csv.writer(f).writerow(METADATA_FIELDS)
        return done_set

    # Chỉ giữ những file thực sự còn tồn tại trên đĩa
    valid_rows = []
    try:
        with open(meta_csv, "r", encoding="utf-8", errors="ignore") as f:
            reader = csv.reader(f)
            next(reader, None)
            for row in reader:
                if row and len(row) >= len(METADATA_FIELDS):
                    fn = row[-1].strip()
                    if is_file_complete(os.path.join(data_dir, fn)):
                        done_set.add(fn)
                        valid_rows.append(row)

        with open(meta_csv, "w", encoding="utf-8", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(METADATA_FIELDS)
            writer.writerows(valid_rows)
    except Exception:
        pass

    return done_set


def main():
    parser = argparse.ArgumentParser(description="Tải và xử lý rút gọn dữ liệu OAS")
    parser.add_argument("--sh", default=r"D:\OAS\bulk_download_human_unpaired.sh", help="Đường dẫn file shell tải dữ liệu")
    parser.add_argument("--data-dir", default=r"D:\OAS\human_unpaired", help="Thư mục lưu file kết quả")
    parser.add_argument("--meta-csv", default=r"D:\OAS\OAS_metadata.csv", help="Đường dẫn lưu metadata CSV")
    parser.add_argument("--threads", type=int, default=8, help="Số luồng tải song song (mặc định 8)")
    parser.add_argument("--fresh", action="store_true", help="Xóa sạch làm lại từ đầu")
    args = parser.parse_args()

    data_dir = args.data_dir
    meta_csv = args.meta_csv
    sh_path = args.sh if os.path.exists(args.sh) else "bulk_download_human_unpaired.sh"

    # Nếu chọn làm mới từ đầu
    if args.fresh:
        print("Đang xóa dữ liệu cũ để chạy lại từ đầu...")
        if os.path.exists(data_dir):
            for f in os.listdir(data_dir):
                try:
                    os.remove(os.path.join(data_dir, f))
                except Exception:
                    pass
        if os.path.exists(meta_csv):
            try:
                os.remove(meta_csv)
            except Exception:
                pass
        print("Đã dọn dẹp xong.\n")

    os.makedirs(data_dir, exist_ok=True)

    # Đọc danh sách metadata và danh sách URL
    existing_meta = load_metadata_set(meta_csv, data_dir)
    all_targets = get_urls_from_script(sh_path)
    total = len(all_targets)

    if total == 0:
        print(f"Không tìm thấy URL nào trong file {sh_path}")
        return

    # Lọc danh sách file còn thiếu
    pending = [item for item in all_targets if not is_file_complete(os.path.join(data_dir, item[0]))]
    already_done = total - len(pending)

    print(f"Tổng số file cần tải:  {total:,}")
    print(f"Đã hoàn thành trước:   {already_done:,}")
    print(f"Cần tải và xử lý tiếp: {len(pending):,}")
    print(f"Số luồng chạy:         {args.threads}")
    print("-" * 75)

    if not pending:
        print("Tất cả các file đã được tải và xử lý xong!")
        return

    t_start = time.time()
    done_count = 0
    total_raw = 0
    total_trimmed = 0

    executor = ThreadPoolExecutor(max_workers=args.threads)
    try:
        futures = {executor.submit(process_url, item, data_dir, meta_csv, existing_meta): item for item in pending}
        for idx, future in enumerate(as_completed(futures), 1):
            fn, status, r_sz, t_sz = future.result()
            elapsed = time.time() - t_start
            speed = idx / elapsed if elapsed > 0 else 0
            eta = format_time((len(pending) - idx) / speed) if speed > 0 else "--:--"

            if status == "OK":
                done_count += 1
                total_raw += r_sz
                total_trimmed += t_sz
                current_total = already_done + done_count
                pct = (current_total / total) * 100
                r_mb = r_sz / (1024 * 1024)
                t_mb = t_sz / (1024 * 1024)
                print(f"[{current_total:5d}/{total:,}] ({pct:5.1f}%) {fn[:32]:32s} | {r_mb:5.1f}MB -> {t_mb:4.1f}MB | ETA: {eta}")
            elif status != "SKIP":
                print(f"Lỗi xử lý file {fn}: {status}")

    except KeyboardInterrupt:
        print("\nĐã nhấn Ctrl+C, đang dừng chương trình...")
        executor.shutdown(wait=False, cancel_futures=True)
        sys.exit(0)
    finally:
        executor.shutdown(wait=True)

    print("-" * 75)
    print(f"Hoàn thành! Thời gian: {format_time(time.time() - t_start)}")
    print(f"Số file mới xử lý: {done_count:,}")
    print(f"Dung lượng raw ước tính: {total_raw / (1024**3):.2f} GB")
    print(f"Dung lượng sau khi nén:  {total_trimmed / (1024**3):.2f} GB")


if __name__ == "__main__":
    main()
