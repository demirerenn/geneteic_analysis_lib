# genetik_analiz.py

import bisect
import collections
import string
import random
import itertools  # YENİ EKLENDİ


# ==============================================================================
# BÖLÜM 1: DOSYA OKUMA VE DNA DİZİSİ İŞLEMLERİ
# ==============================================================================

def readGenome(filename):
    """
    Bir FASTA dosyasını okur ve genom dizisini tek bir string olarak döndürür.
    Header satırlarını ('>' ile başlayanlar) yok sayar.

    Args:
        filename (str): Okunacak FASTA dosyasının adı.

    Returns:
        str: Genom dizisi.
    """
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                genome += line.rstrip()
    return genome


def readFastq(filename):
    """
    Bir FASTQ dosyasını okur ve dizileri (sequences) ile kalite skorlarını (qualities)
    iki ayrı liste olarak döndürür.

    Args:
        filename (str): Okunacak FASTQ dosyasının adı.

    Returns:
        tuple: (sequences, qualities) içeren bir tuple.
    """
    sequences = []
    qualities = []
    with open(filename) as fh:
        while True:
            fh.readline()  # name line
            seq = fh.readline().rstrip()
            fh.readline()  # placeholder line
            qual = fh.readline().rstrip()
            if len(seq) == 0:
                break
            sequences.append(seq)
            qualities.append(qual)
    return sequences, qualities


def reverseComplement(s):
    """
    Bir DNA dizisinin ters tamamlayıcısını (reverse complement) döndürür.

    Args:
        s (str): DNA dizisi.

    Returns:
        str: Dizinin ters tamamlayıcısı.
    """
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    t = ''
    for base in s:
        t = complement[base] + t
    return t


def phred33ToQ(qual):
    """
    Phred+33 formatındaki bir kalite karakterini tamsayıya çevirir.

    Args:
        qual (str): Tek bir kalite karakteri.

    Returns:
        int: Phred kalite skoru.
    """
    return ord(qual) - 33


# ==============================================================================
# BÖLÜM 2: BASİT EŞLEŞTİRME VE ANALİZ FONKSİYONLARI
# ==============================================================================

def naive(p, t):
    """
    Naif (brute-force) tam eşleştirme algoritması.

    Args:
        p (str): Aranacak kalıp (pattern).
        t (str): Aranacak metin (text/genome).

    Returns:
        list: Eşleşmelerin başlangıç indislerini içeren liste.
    """
    occurrences = []
    for i in range(len(t) - len(p) + 1):
        match = True
        for j in range(len(p)):
            if t[i + j] != p[j]:
                match = False
                break
        if match:
            occurrences.append(i)
    return occurrences


def findGCByPos(reads):
    """
    Read'lerin her pozisyonundaki GC oranını hesaplar.

    Args:
        reads (list): DNA read dizilerinin listesi.

    Returns:
        list: Her pozisyon için hesaplanmış GC oranını içeren liste.
    """
    gc = [0] * 100
    totals = [0] * 100
    for read in reads:
        for i in range(len(read)):
            if read[i] == 'C' or read[i] == 'G':
                gc[i] += 1
            totals[i] += 1
    for i in range(len(gc)):
        if totals[i] > 0:
            gc[i] /= float(totals[i])
    return gc


def generateReads(genome, numReads, readLen):
    """
    Verilen bir genomdan rastgele pozisyonlardan read'ler üretir.

    Args:
        genome (str): Genom dizisi.
        numReads (int): Üretilecek read sayısı.
        readLen (int): Her bir read'in uzunluğu.

    Returns:
        list: Üretilen read'lerin listesi.
    """
    reads = []
    for _ in range(numReads):
        start = random.randint(0, len(genome) - readLen)
        reads.append(genome[start: start + readLen])
    return reads


def editDistance(x, y):
    """
    Dinamik programlama kullanarak iki string arasındaki Levenshtein (edit)
    mesafesini hesaplar.

    Args:
        x (str): Birinci string.
        y (str): İkinci string.

    Returns:
        int: İki string arasındaki edit mesafesi.
    """
    D = []
    for i in range(len(x) + 1):
        D.append([0] * (len(y) + 1))

    for i in range(len(x) + 1):
        D[i][0] = i
    for j in range(len(y) + 1):
        D[0][j] = j

    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            distHor = D[i][j - 1] + 1
            distVer = D[i - 1][j] + 1
            if x[i - 1] == y[j - 1]:
                distDiag = D[i - 1][j - 1]
            else:
                distDiag = D[i - 1][j - 1] + 1
            D[i][j] = min(distHor, distVer, distDiag)

    return D[-1][-1]


# ==============================================================================
# BÖLÜM 3: GELİŞMİŞ EŞLEŞTİRME VE HİZALAMA ALGORİTMALARI
# ==============================================================================

# ------------------------------------------------------------------------------
# Boyer-Moore Algoritması için Yardımcı Fonksiyonlar ve Sınıf
# ------------------------------------------------------------------------------

def z_array(s):
    assert len(s) > 1
    z = [len(s)] + [0] * (len(s) - 1)
    for i in range(1, len(s)):
        if s[i] == s[i - 1]:
            z[1] += 1
        else:
            break
    r, l = 0, 0
    if z[1] > 0:
        r, l = z[1], 1
    for k in range(2, len(s)):
        assert z[k] == 0
        if k > r:
            for i in range(k, len(s)):
                if s[i] == s[i - k]:
                    z[k] += 1
                else:
                    break
            r, l = k + z[k] - 1, k
        else:
            nbeta = r - k + 1
            zkp = z[k - l]
            if nbeta > zkp:
                z[k] = zkp
            else:
                nmatch = 0
                for i in range(r + 1, len(s)):
                    if s[i] == s[i - k]:
                        nmatch += 1
                    else:
                        break
                l, r = k, r + nmatch
                z[k] = r - k + 1
    return z


def n_array(s):
    return z_array(s[::-1])[::-1]


def big_l_prime_array(p, n):
    lp = [0] * len(p)
    for j in range(len(p) - 1):
        i = len(p) - n[j]
        if i < len(p):
            lp[i] = j + 1
    return lp


def big_l_array(p, lp):
    l = [0] * len(p)
    l[1] = lp[1]
    for i in range(2, len(p)):
        l[i] = max(l[i - 1], lp[i])
    return l


def small_l_prime_array(n):
    small_lp = [0] * len(n)
    for i in range(len(n)):
        if n[i] == i + 1:
            small_lp[len(n) - i - 1] = i + 1
    for i in range(len(n) - 2, -1, -1):
        if small_lp[i] == 0:
            small_lp[i] = small_lp[i + 1]
    return small_lp


def good_suffix_table(p):
    n = n_array(p)
    lp = big_l_prime_array(p, n)
    return lp, big_l_array(p, lp), small_l_prime_array(n)


def dense_bad_char_tab(p, amap):
    tab = []
    nxt = [0] * len(amap)
    for i in range(0, len(p)):
        c = p[i]
        assert c in amap
        tab.append(nxt[:])
        nxt[amap[c]] = i + 1
    return tab


class BoyerMoore(object):
    """ Boyer-Moore algoritmasını uygular. """

    def __init__(self, p, alphabet='ACGT'):
        self.p = p
        self.alphabet = alphabet
        self.amap = {self.alphabet[i]: i for i in range(len(self.alphabet))}
        self.bad_char = dense_bad_char_tab(p, self.amap)
        _, self.big_l, self.small_l_prime = good_suffix_table(p)

    def bad_character_rule(self, i, c):
        ci = self.amap[c]
        return i - (self.bad_char[i][ci] - 1)

    def good_suffix_rule(self, i):
        length = len(self.big_l)
        if i == length - 1:
            return 0
        i += 1
        if self.big_l[i] > 0:
            return length - self.big_l[i]
        return length - self.small_l_prime[i]

    def match_skip(self):
        return len(self.small_l_prime) - self.small_l_prime[1]


def boyer_moore(p, t):
    """
    Boyer-Moore algoritmasını kullanarak tam eşleştirme yapar.

    Args:
        p (str): Aranacak kalıp (pattern).
        t (str): Aranacak metin (text/genome).

    Returns:
        list: Eşleşmelerin başlangıç indislerini içeren liste.
    """
    p_bm = BoyerMoore(p, alphabet='ACGTN')  # 'N' de ekleyebiliriz
    i = 0
    occurrences = []
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        for j in range(len(p) - 1, -1, -1):
            if p[j] != t[i + j]:
                skip_bc = p_bm.bad_character_rule(j, t[i + j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences


# ------------------------------------------------------------------------------
# Substring Index (k-mer Index) ile Eşleştirme
# ------------------------------------------------------------------------------

class Index(object):
    """ Belirtilen k uzunluğundaki tüm alt diziler için bir indeks oluşturur. """

    def __init__(self, t, k):
        self.k = k
        self.index = []
        for i in range(len(t) - k + 1):
            self.index.append((t[i:i + k], i))
        self.index.sort()

    def query(self, p):
        """ Kalıbın ilk k-mer'i için indeksteki eşleşmeleri döndürür. """
        kmer = p[:self.k]
        i = bisect.bisect_left(self.index, (kmer, -1))
        hits = []
        while i < len(self.index):
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits


def queryIndex(p, t, index):
    """
    Bir k-mer indeksi kullanarak tam eşleştirme yapar.

    Args:
        p (str): Aranacak kalıp (pattern).
        t (str): Aranacak metin (text/genome).
        index (Index): Önceden oluşturulmuş Index nesnesi.

    Returns:
        list: Eşleşmelerin başlangıç indislerini içeren liste.
    """
    k = index.k
    offsets = []
    for i in index.query(p):
        if p[k:] == t[i + k:i + len(p)]:
            offsets.append(i)
    return offsets


# ------------------------------------------------------------------------------
# Yaklaşık Eşleştirme (Approximate Matching)
# ------------------------------------------------------------------------------

def approximate_match(p, t, n):
    """
    Pigeonhole prensibini kullanarak en fazla 'n' uyuşmazlığa (mismatch)
    izin veren yaklaşık eşleştirme yapar.

    Args:
        p (str): Aranacak kalıp (pattern).
        t (str): Aranacak metin (text/genome).
        n (int): İzin verilen maksimum uyuşmazlık sayısı.

    Returns:
        list: Yaklaşık eşleşmelerin başlangıç indislerini içeren liste.
    """
    segment_length = int(round(len(p) / (n + 1)))
    all_matches = set()
    for i in range(n + 1):
        start = i * segment_length
        end = min((i + 1) * segment_length, len(p))
        p_segment = p[start:end]

        p_bm = BoyerMoore(p_segment, alphabet='ACGTN')
        matches = boyer_moore(p_segment, t)

        for m in matches:
            if m < start or m - start + len(p) > len(t):
                continue
            mismatches = 0
            for j in range(0, start):
                if not p[j] == t[m - start + j]:
                    mismatches += 1
            for j in range(end, len(p)):
                if not p[j] == t[m - start + j]:
                    mismatches += 1

            if mismatches <= n:
                all_matches.add(m - start)

    return list(all_matches)


# ------------------------------------------------------------------------------
# Global Alignment (Needleman-Wunsch)
# ------------------------------------------------------------------------------

# Skorlama matrisi: A, C, G, T ve boşluk (-) için uyum/uyumsuzluk/boşluk cezaları
alphabet = ['A', 'C', 'G', 'T']
score = [[0, 4, 2, 4, 8],  # A ile A,C,G,T,- skorları
         [4, 0, 4, 2, 8],  # C ile A,C,G,T,- skorları
         [2, 4, 0, 4, 8],  # G ile A,C,G,T,- skorları
         [4, 2, 4, 0, 8],  # T ile A,C,G,T,- skorları
         [8, 8, 8, 8, 8]]  # - (boşluk) ile A,C,G,T,- skorları


def globalAlignment(x, y):
    """
    Needleman-Wunsch algoritmasını kullanarak iki dizi arasındaki en iyi
    global hizalama skorunu hesaplar. Minimum skor hedeflenir.

    Args:
        x (str): Birinci dizi.
        y (str): İkinci dizi.

    Returns:
        int: Optimal global hizalama skoru.
    """
    # Uzaklık matrisini oluştur
    D = []
    for i in range(len(x) + 1):
        D.append([0] * (len(y) + 1))

    # İlk sütunu başlat
    for i in range(1, len(x) + 1):
        D[i][0] = D[i - 1][0] + score[alphabet.index(x[i - 1])][-1]

    # İlk satırı başlat
    for j in range(1, len(y) + 1):
        D[0][j] = D[0][j - 1] + score[-1][alphabet.index(y[j - 1])]

    # Matrisin geri kalanını doldur
    for i in range(1, len(x) + 1):
        for j in range(1, len(y) + 1):
            distHor = D[i][j - 1] + score[-1][alphabet.index(y[j - 1])]
            distVer = D[i - 1][j] + score[alphabet.index(x[i - 1])][-1]
            distDiag = D[i - 1][j - 1] + score[alphabet.index(x[i - 1])][alphabet.index(y[j - 1])]
            D[i][j] = min(distHor, distVer, distDiag)

    # Sağ alt köşedeki değeri döndür
    return D[-1][-1]


# ==============================================================================
# BÖLÜM 4: DİZİ BİRLEŞTİRME (ASSEMBLY) İÇİN OVERLAP FONKSİYONLARI - YENİ
# ==============================================================================

def overlap(a, b, min_length=3):
    """
    'a' dizisinin 'b' dizisinin başına uyan en uzun son parçasının (suffix)
    uzunluğunu döndürür. Bu parçanın uzunluğu en az 'min_length' olmalıdır.
    Eğer böyle bir overlap yoksa 0 döndürür.

    Args:
        a (str): Birinci dizi.
        b (str): İkinci dizi.
        min_length (int): Minimum overlap uzunluğu.

    Returns:
        int: En uzun overlap'in uzunluğu veya 0.
    """
    start = 0
    while True:
        # b'nin ilk min_length parçasını a içinde ara
        start = a.find(b[:min_length], start)
        if start == -1:  # Daha fazla eşleşme yoksa
            return 0
        # Eşleşme bulundu; tam suffix/prefix uyumunu kontrol et
        if b.startswith(a[start:]):
            return len(a) - start
        start += 1  # Bir sonraki pozisyondan aramaya devam et


def naive_overlap_map(reads, k):
    """
    Verilen read'ler arasındaki tüm olası overlap'leri bulur.

    Args:
        reads (list): DNA read'lerinin listesi.
        k (int): Minimum overlap uzunluğu.

    Returns:
        dict: (read1, read2): overlap_uzunlugu formatında bir sözlük.
    """
    olaps = {}
    for a, b in itertools.permutations(reads, 2):
        olen = overlap(a, b, min_length=k)
        if olen > 0:
            olaps[(a, b)] = olen
    return olaps