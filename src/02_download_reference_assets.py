"""Download large public reference assets with resumable HTTP range requests."""

from __future__ import annotations

import argparse
import concurrent.futures
import time
import urllib.request
from pathlib import Path


def remote_size(url: str) -> int:
    request = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(request, timeout=60) as response:
        length = response.headers.get("Content-Length")
    if length is None:
        raise RuntimeError(f"Server did not report Content-Length: {url}")
    return int(length)


def download_range(url: str, output: Path, start: int, end: int, retries: int) -> int:
    expected = end - start + 1
    for attempt in range(1, retries + 1):
        try:
            request = urllib.request.Request(url, headers={"Range": f"bytes={start}-{end}"})
            with urllib.request.urlopen(request, timeout=120) as response:
                if response.status != 206:
                    raise RuntimeError(f"Range request returned HTTP {response.status}")
                data = response.read()
            if len(data) != expected:
                raise RuntimeError(f"Expected {expected} bytes, received {len(data)}")
            with output.open("r+b", buffering=0) as handle:
                handle.seek(start)
                handle.write(data)
            return expected
        except Exception:
            if attempt == retries:
                raise
            time.sleep(min(2**attempt, 10))
    raise AssertionError("unreachable")


def download(url: str, output: Path, workers: int, retries: int) -> None:
    size = remote_size(url)
    output.parent.mkdir(parents=True, exist_ok=True)
    marker = output.with_name(output.name + ".complete")
    progress = output.with_name(output.name + ".ranges")
    if marker.exists() and output.exists() and output.stat().st_size == size:
        print(f"SKIP complete {output} ({size} bytes)", flush=True)
        return

    if not output.exists() or output.stat().st_size != size or not progress.exists():
        with output.open("wb") as handle:
            handle.truncate(size)
        progress.write_text("", encoding="ascii")

    chunk = max(4 * 1024 * 1024, (size + workers - 1) // workers)
    ranges = [(start, min(start + chunk - 1, size - 1)) for start in range(0, size, chunk)]
    done = set(progress.read_text(encoding="ascii").splitlines())
    completed = sum(end - start + 1 for start, end in ranges if f"{start}-{end}" in done)
    pending = [(start, end) for start, end in ranges if f"{start}-{end}" not in done]
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(download_range, url, output, start, end, retries): (start, end)
            for start, end in pending
        }
        for future in concurrent.futures.as_completed(futures):
            start, end = futures[future]
            completed += future.result()
            with progress.open("a", encoding="ascii") as handle:
                handle.write(f"{start}-{end}\n")
            print(f"{output.name}: {completed / size:.1%}", flush=True)

    if output.stat().st_size != size:
        raise RuntimeError(f"Size mismatch after download: {output}")
    marker.write_text(f"{url}\n{size}\n", encoding="utf-8")
    progress.unlink(missing_ok=True)
    print(f"DONE {output} ({size} bytes)", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("url")
    parser.add_argument("output", type=Path)
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--retries", type=int, default=5)
    args = parser.parse_args()
    if not 1 <= args.workers <= 64:
        parser.error("--workers must be between 1 and 64")
    download(args.url, args.output, args.workers, args.retries)


if __name__ == "__main__":
    main()
