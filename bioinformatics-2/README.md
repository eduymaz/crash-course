# Snakemake Workflow Tutorial

Bu klasör, Snakemake ile basit bir veri işleme pipeline'ı nasıl oluşturacağınızı anlatan rehber niteliğinde örnekler içerir.

## 📦 Gereksinimler

- Python 3.12+
- Snakemake (kurulum önerisi: `conda install -c bioconda snakemake`)

## 📁 Dizin Yapısı

```
bioinformatics-2/
├── sample/                # Girdi dosyaları
│   ├── sample1.txt
│   └── sample2.txt
├── scripts/               # Yardımcı Python/Komut dosyaları
│   └── uppercase.py       # Örnek Python scripti
├── Snakefile              # Snakemake workflow tanımı
└── README.md              # Bu rehber
```

## ⚙️ Kurulum

```bash
# Conda ile kurulum (tercih edilen yöntem)
conda create -n snakemake_env snakemake python=3.12.8
conda activate snakemake_env
```

## 🛠️ Basit Workflow Örneği

Bu örnek, her bir sample dosyasındaki satır sayısını hesaplayıp `results/` dizinine yazar.

### Snakefile İçeriği
```python
rule all:
    input:
        expand("results/{sample}.count", sample=["sample1","sample2"])

rule count_lines:
    input:
        "sample/{sample}.txt"
    output:
        "results/{sample}.count"
    shell:
        "wc -l {input} > {output}"
```

### Nasıl Çalıştırılır?
```bash
snakemake --cores 1 #cihaz gücünüz ve işlem yükünüze göre 12 24 şeklinde artırılabilir 
```

## 🐍 Python Script Kullanımı

`uppercase.py` scripti, girdi dosyasını büyük harfe çevirip çıktıyı yazar:

```python
import sys

def main(input_file, output_file):
    with open(input_file) as inf, open(output_file, 'w') as outf:
        data = inf.read().upper()
        outf.write(data)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
```

Snakefile'de bu scripti şöyle çağırabilirsiniz:
```python
rule uppercase:
    input:
        "sample/{sample}.txt"
    output:
        "results/{sample}.upper.txt"
    script:
        "scripts/uppercase.py"
```

## 📊 Çıktı Klasörü

`results/` altında şu dosyalar oluşacaktır:

- `sample1.count`, `sample2.count` 
- `sample1.upper.txt`, `sample2.upper.txt`

## 🔬 Gelişmiş Workflow: Mapping ve Variant Calling

Aşağıdaki örnek, temel satır sayısı ve uppercase kurallarına ek olarak:
1. Reads'in referans genom ile hizalanması (BWA)
2. Hizalanan BAM dosyasının sıralanması ve indekslenmesi (samtools)
3. Varyant çağrımı (bcftools)

aşağıdaki gibi tanımlanır:

```python
rule all:
    input:
        # Basit çıktılar
        expand("results/{sample}.count", sample=["sample1","sample2"]),
        expand("results/{sample}.upper.txt", sample=["sample1","sample2"]),
        # Gelişmiş çıktılar
        expand("alignments/{sample}.sorted.bam", sample=["sample1","sample2"]),
        expand("variants/{sample}.vcf", sample=["sample1","sample2"])

# 1) Reads hizalama
rule map_reads:
    input:
        reads="sample/{sample}.txt",
        ref="data/ref.fa"
    output:
        bam="alignments/{sample}.bam"
    shell:
        "bwa mem {input.ref} {input.reads} | samtools view -Sb - > {output.bam}"

# 2) BAM sıralama
rule sort_bam:
    input:
        bam="alignments/{sample}.bam"
    output:
        sorted_bam="alignments/{sample}.sorted.bam"
    shell:
        "samtools sort {input.bam} -o {output.sorted_bam}"

# 3) BAM indeksleme
rule index_bam:
    input:
        sorted_bam="alignments/{sample}.sorted.bam"
    output:
        bai="alignments/{sample}.sorted.bam.bai"
    shell:
        "samtools index {input.sorted_bam}"

# 4) Varyant çağrımı
rule call_variants:
    input:
        bam="alignments/{sample}.sorted.bam",
        bai="alignments/{sample}.sorted.bam.bai",
        ref="data/ref.fa"
    output:
        vcf="variants/{sample}.vcf"
    shell:
        "bcftools mpileup -f {input.ref} {input.bam} | bcftools call -mv -Ov -o {output.vcf}"
```

**Not:** Bu adımları kullanabilmek için `data/` klasörüne `ref.fa` (referans genom FASTA) dosyası koymalısınız. Snakemake, her kuralı otomatik olarak doğru sırayla çalıştırarak bağımlılık grafiğini (DAG) oluşturur ve paralel çalışmayı yönetir.

---

Bu adımları takip ederek karmaşık işlem adımlarını birbirine bağlayan, ölçeklenebilir bir Snakemake pipeline'ı geliştirebilirsiniz.