# NF-Core Crash Course: Sarek Pipeline v3.5.1

Bu rehber, NF-Core altyapısını kullanarak **Sarek** exome/whole-genome sekans analizi pipeline'ını nasıl başlatıp çalıştıracağınızı adım adım gösterir.

## 🛠️ 1) Gereksinimler

- Java 8+ (Nextflow için)
- Nextflow (>=20.07.1)
- nf-core/tools (CLI yardımcı aracı)
- Docker veya Singularity (container çalıştırmak için)

```bash
java -version

curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/

pip install nf-core
```

## 📥 2) Pipeline Oluşturma/İndirme

NF-Core ile doğrudan `nf-core/sarek` pipeline'ını çalıştırmak için manual download gerekmez, `nextflow run` komutu her şeyi çeker:

```bash
tcd bioinformatics-1

nextflow run nf-core/sarek \
    -r 3.5.1 \
    --input sample/reads/ \
    --genome GRCh38 \
    --outdir sarek_results/ \
    --max_cpus 8 \
    --max_memory '16.GB' \
    -profile docker
```

- `-r 3.5.1`: Versiyonu sabitler (dokümantasyondaki örnek)
- `--input`: FASTQ dosyalarının bulunduğu klasör (`sample/reads/*.fastq.gz`)
- `--genome`: Referans (örn. GRCh37, hg38)
- `--outdir`: Çıktı klasörü
- `-profile`: `docker`, `singularity` veya `conda`

## 📂 2b) YAML Konfigürasyon ile Çalıştırma

Eğer `sarek_run.yaml` gibi bir config dosyası hazırladıysanız, parametreleri komut satırında vermek yerine bu dosyayı kullanabilirsiniz:

```bash
nextflow run nf-core/sarek \
    -r 3.5.1 \
    -c sarek_run.yaml
```

Bu yöntemle `--input`, `--genome`, `--outdir`, `--max_cpus` ve örneğin `--samplesheet` gibi tüm parametreleri `sarek_run.yaml` içinde bir arada saklayarak çalışma sürecini daha okunabilir ve yeniden üretilebilir hale getirirsiniz.

## ⚙️ 3) Parametre Özelleştirme

`nextflow run` komutuna ek olarak birçok parametre ekleyebilirsiniz:
```bash
--bwa_index <path/to/index>     # Özel BWA dizin kastlama
--vcf_filter <expression>       # VCF varyant filtresi
--skip_qc                       # QC adımını atla
--publish_dir <dir>:format("fq.gz")
```
Ayrıntılı parametre listesi için:
```bash
nextflow run nf-core/sarek -r 3.5.1 --help
```

## 📊 4) Çıktılar

- `sarek_results/fastqc/`: Raw QC raporları
- `sarek_results/bwa/`: Hizalanmış BAM dosyaları
- `sarek_results/gatk/`: VCF varyant dosyaları ve filtrelenmiş sonuçlar
- Log, MultiQC raporları ve metrikler

```bash
grep "#CHROM" sarek_results/gatk/*.vcf | head -n 20
```

## 🧠 5) Crash Course İpuçları

1. İlk denemede `-profile docker` daha temizdir; `conda` profilini kullanmak için `-profile conda` ekleyin.  
2. Pipeline'ı ara adımlarla debug etmek için `-with-report report.html` ve `-with-trace trace.txt` kullanın.  
3. Paralelleştirmeyi `--max_cpus` ile kontrol ederek makinada aşırı yüklenmeyi önleyin.  
4. Kanal (channel) yapılarını anlamak ve özel kural eklemek için Nextflow DSL2 belgelerine bakın.

---

*Happy sequencing with NF-Core!* 