import os

def read_stream(path: str) -> str:
    with open(path, "r", encoding="utf-8") as f:
        return "".join(line.strip() for line in f)
    

def z_algorithm(s: str):
    n = len(s)
    Z = [0] * n
    l = r = 0
    for i in range(1, n):
        if i <= r:
            Z[i] = min(r - i + 1, Z[i - 1])
        while i + Z[i] < n and s[Z[i]] == s[i + Z[i]]:
            Z[i] += 1
        if i + Z[i] - 1 > r:
            l, r = i, i + Z[i] - 1
    return Z


def pattern_z(text: str, pattern: str):
    if not pattern:
        return True, 1
    
    concat = pattern + "$" + text
    Z = z_algorithm(concat)
    pat_len = len(pattern)
    first_position = None

    for i in range(pat_len + 1, len(concat)):
        if Z[i] == pat_len:
            pos0 = i - (pat_len + 1)
            first_position = pos0 + 1
            break

    if first_position is None:
        return False, None
    else:
        return True, first_position
    

def longest_palindrome_pos(s: str):
    n = len(s)
    if n == 0:
        return 1, 0
    
    best_len = 1
    best_start = 0

    def expand(left: int, right: int):
        nonlocal best_len, best_start
        while left >= 0 and right < n and s[left] == s[right]:
            curr_len = right - left + 1
            if curr_len > best_len:
                best_len = curr_len
                best_start = left
            left -= 1
            right += 1

    for center in range(n):
        expand(center, center)
        expand(center, center + 1)

    start_1 = best_start + 1
    end_1 = best_start + best_len
    return start_1, end_1


def lcs_pos(s: str, t: str):
    n, m = len(s), len(t)
    if n == 0 or m == 0:
        return 1, 0
    
    prev = [0] * (m + 1)
    curr = [0] * (m + 1)

    best_len = 0
    best_end_i = 0

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if s[i - 1] == t[j - 1]:
                curr[j] = prev[j - 1] + 1
                if curr[j] > best_len:
                    best_len = curr[j]
                    best_end_i = i
            else:
                curr[j] = 0

        prev, curr = curr, prev

    if best_len == 0:
        return 1, 0
    
    start_1 = best_end_i - best_len + 1
    end_1 = best_end_i
    return start_1, end_1


def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    data_dir = os.path.join(base_dir, "data")

    transmission1 = read_stream(os.path.join(data_dir, "transmission1.txt"))
    transmission2 = read_stream(os.path.join(data_dir, "transmission2.txt"))
    mcode1 = read_stream(os.path.join(data_dir, "mcode1.txt"))
    mcode2 = read_stream(os.path.join(data_dir, "mcode2.txt"))
    mcode3 = read_stream(os.path.join(data_dir, "mcode3.txt"))

    for text in (transmission1, transmission2):
        for mcode in (mcode1, mcode2, mcode3):
            found, pos = pattern_z(text, mcode)
            if found:
                print(f"true {pos}")
            else:
                print("false")

    start1, end1 = longest_palindrome_pos(transmission1)
    start2, end2 = longest_palindrome_pos(transmission2)
    print(f"{start1} {end1}")
    print(f"{start2} {end2}")

    start_lcs, end_lcs = lcs_pos(transmission1, transmission2)
    print(f"{start_lcs} {end_lcs}")

if __name__ == "__main__":
    main()