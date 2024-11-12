# Changelog

## [1.2.0](https://github.com/3d-omics/hg_genotype/compare/v1.1.4...v1.2.0) (2024-11-12)


### Features

* demand snakemake 8 ([10030be](https://github.com/3d-omics/hg_genotype/commit/10030be50554e148e25b647ff2fedbe7f5b9fdff))


### Bug Fixes

* add all the gzis and fais as necessary ([a68f5af](https://github.com/3d-omics/hg_genotype/commit/a68f5af678a74cc88b517ecfbc1dc52041a479c9))
* add gzi and fai ([3e647f3](https://github.com/3d-omics/hg_genotype/commit/3e647f32f3336653ad2b934602419469c2b64f87))
* add misssing file ([b00362a](https://github.com/3d-omics/hg_genotype/commit/b00362a329722f8b7ef1d91a5db2bcc3d861fbef))
* ask for fai and gzi ([9aedadf](https://github.com/3d-omics/hg_genotype/commit/9aedadf9cad9c39d1768c78d716518410d644820))
* don't get stats from  snakehelpers ([38b7b5c](https://github.com/3d-omics/hg_genotype/commit/38b7b5cb6fed4cb9a84438505def0ec098175593))
* execute multiqc swaps ([3aeb3a6](https://github.com/3d-omics/hg_genotype/commit/3aeb3a6b1bd902aaf6615e3799dde1fb46ca5356))
* typo ([a5ce9a5](https://github.com/3d-omics/hg_genotype/commit/a5ce9a5f9904ab437586127c008db1bb7901698f))
* wrong wildcard names ([53b9dcc](https://github.com/3d-omics/hg_genotype/commit/53b9dcc1827a9a838b64541ec077e0e34dd0155f))


### Performance Improvements

* group align per sample ([c7392f7](https://github.com/3d-omics/hg_genotype/commit/c7392f79e3ea69a4c494f441ad348d5347f2b4e7))
* omit caching because of version in bwamem2 and reference ([15b3a9d](https://github.com/3d-omics/hg_genotype/commit/15b3a9d8613a49579a954e5267e99217e22936ef))

## [1.1.4](https://github.com/3d-omics/hg_genotype/compare/v1.1.3...v1.1.4) (2024-08-06)


### Bug Fixes

* upload schema.svg ([675f4de](https://github.com/3d-omics/hg_genotype/commit/675f4dec7d60fc4ee4245a8231f36a4fa258bed6))

## [1.1.3](https://github.com/3d-omics/hg_genotype/compare/v1.1.2...v1.1.3) (2024-08-06)


### Bug Fixes

* correct rule names ([d002656](https://github.com/3d-omics/hg_genotype/commit/d0026563ddca0a50ef966c37c5eb938f1ff60449))
* disable groups ([2461f2e](https://github.com/3d-omics/hg_genotype/commit/2461f2e1f2b753265776279a499d0c1cec2e44dd))
* path to string ([8de31b0](https://github.com/3d-omics/hg_genotype/commit/8de31b025663ca43827e5ca05e750c8a4e1811d5))


### Performance Improvements

* update profile with better values ([9c9e80e](https://github.com/3d-omics/hg_genotype/commit/9c9e80e74c652bdcb4f75aeeb415ac63384b991b))

## [1.1.2](https://github.com/3d-omics/hg_genotype/compare/v1.1.1...v1.1.2) (2024-06-06)


### Bug Fixes

* convert param from path to string ([1cb8ca3](https://github.com/3d-omics/hg_genotype/commit/1cb8ca39bf10944900d4d97fea4914bab86238f1))
* match resources between piped rules ([19b5a41](https://github.com/3d-omics/hg_genotype/commit/19b5a41cf66f66b8c179619697f5c8e00742df47))

## [1.1.1](https://github.com/3d-omics/hg_genotype/compare/v1.1.0...v1.1.1) (2024-06-06)


### Bug Fixes

* infer properly sex chromosomes when only one is in the reference ([e89e48c](https://github.com/3d-omics/hg_genotype/commit/e89e48c110c9f6f5d5809105a7d1695ac8656194))
* put all rule on top of cache ([d404cb6](https://github.com/3d-omics/hg_genotype/commit/d404cb65b5239022db2c4a17be18b7f37d379bc5))

## [1.1.0](https://github.com/3d-omics/hg_genotype/compare/v1.0.0...v1.1.0) (2024-05-13)


### Features

* add bwa index to cache ([33d45a4](https://github.com/3d-omics/hg_genotype/commit/33d45a4e5f3ff2dee058f5f521a874c6373afad1))
* add cache rule ([f64131f](https://github.com/3d-omics/hg_genotype/commit/f64131f5125790a7c03e94f448d05f8ecde795c3))
* add cache to reference rules ([2ebb95d](https://github.com/3d-omics/hg_genotype/commit/2ebb95d7c4d9abf7ee8f6f02bbd508f6a3a132eb))


### Bug Fixes

* correct input for metarule align__index ([37fbc97](https://github.com/3d-omics/hg_genotype/commit/37fbc972ddcc9fe74d79b2e49280ec0fdc731ddb))

## 1.0.0 (2024-05-10)


### Features

* add align and genotype groups ([9510d65](https://github.com/3d-omics/hg_genotype/commit/9510d654ee16204bcca55e6253ac717f0cea3b2c))
* add bedtools ([9df30a5](https://github.com/3d-omics/hg_genotype/commit/9df30a525ec445d3d890fd2aaaa0dbef5f051a8a))
* add groups ([462e389](https://github.com/3d-omics/hg_genotype/commit/462e3897828ad32b5fb6e5eb5ce96509e67b17de))
* add more chromosomes to mock data ([c500757](https://github.com/3d-omics/hg_genotype/commit/c500757595a824d4165eaa0a290d22bce8f155d9))
* add preliminary bcftools calls for gatk ([7db0528](https://github.com/3d-omics/hg_genotype/commit/7db0528b8e414a9485bc175b2981f58db9c78036))
* add resources to the new rules that seem to work at least in a macro level ([30829d2](https://github.com/3d-omics/hg_genotype/commit/30829d2fbae25fcfad546ba01c66e53fda452006))
* add somalier execution + reports ([0ea0640](https://github.com/3d-omics/hg_genotype/commit/0ea06408dc480e85e0f88391ca78952bcca3329f))
* add swaps checkpoint ([80e7758](https://github.com/3d-omics/hg_genotype/commit/80e7758a30a0da82d557a25b9c91a8a23a3d71fc))
* Added sommalier steps. Will try to figure it out how it works in mice ([#8](https://github.com/3d-omics/hg_genotype/issues/8)) ([3ceb0a9](https://github.com/3d-omics/hg_genotype/commit/3ceb0a967461906bb83230bd9ef0eda45a1ce010))
* bowtie2 -&gt; bwa ([005cc35](https://github.com/3d-omics/hg_genotype/commit/005cc35bd3d7283c01b7e7067c4088586b971247))
* convert markduplicates to cram ([a282483](https://github.com/3d-omics/hg_genotype/commit/a28248352272d71078089d819627d42b2b7980fe))
* disable swaps ([f008602](https://github.com/3d-omics/hg_genotype/commit/f008602acf9a0b3226afda5761cef8efe02d7932))
* do vep annotation by sample ([a7dfd4d](https://github.com/3d-omics/hg_genotype/commit/a7dfd4d965f445451babb72420ea3851da296196))
* Document all functions and rules ([#10](https://github.com/3d-omics/hg_genotype/issues/10)) ([265aeb6](https://github.com/3d-omics/hg_genotype/commit/265aeb6b4f4316972f8151302760397827fee2b5))
* Fix versions ([#4](https://github.com/3d-omics/hg_genotype/issues/4)) ([48689dc](https://github.com/3d-omics/hg_genotype/commit/48689dce58f0a1ea20fdfac1056cf204bb256373))
* gtcheck, fastp reports, cache misses ([#18](https://github.com/3d-omics/hg_genotype/issues/18)) ([ca99476](https://github.com/3d-omics/hg_genotype/commit/ca994764cd21b970d317de9b4a5b085284041018))
* i don't remember ([f82ed90](https://github.com/3d-omics/hg_genotype/commit/f82ed906de8760b0a3a5d5bc631b18a87c1827e1))
* little turds ([#13](https://github.com/3d-omics/hg_genotype/issues/13)) ([5920444](https://github.com/3d-omics/hg_genotype/commit/5920444fef43f94e7a6d5e2acfacf32a2fc3aa8f))
* Little turds part 2 ([#14](https://github.com/3d-omics/hg_genotype/issues/14)) ([cac2932](https://github.com/3d-omics/hg_genotype/commit/cac293261500f64f7521be4ac68308218834f5f5))
* merge align bams into a single cram ([2d8d4bc](https://github.com/3d-omics/hg_genotype/commit/2d8d4bc06d101480661ed77061048a65394c446c))
* merge align reports by sample-library, rather than sample-library-chromosome ([0f1a550](https://github.com/3d-omics/hg_genotype/commit/0f1a550bc7bf88e3246dd4ea3ad0981d580dc66d))
* merge libraries into sample in markduplicates ([323837b](https://github.com/3d-omics/hg_genotype/commit/323837b34e4943f46b63f152bd0571a493689158))
* merge libraries into sample in reports ([cd985fd](https://github.com/3d-omics/hg_genotype/commit/cd985fd7b0a36c20c87d70a4e49131f322353e27))
* minimum viable product ([#1](https://github.com/3d-omics/hg_genotype/issues/1)) ([407fa5d](https://github.com/3d-omics/hg_genotype/commit/407fa5d7e93271f7e75573286c96fc6cf96a20ec))
* picard and gatk analyses and reports by chromosome ([#7](https://github.com/3d-omics/hg_genotype/issues/7)) ([d443110](https://github.com/3d-omics/hg_genotype/commit/d443110943dd685e8a3c3e2157d1ce10d6986cd9))
* pin bwa version ([0df670a](https://github.com/3d-omics/hg_genotype/commit/0df670ab121c4796829ace70dbe622333ac72686))
* process samples in haplotype caller ([363cbcd](https://github.com/3d-omics/hg_genotype/commit/363cbcd6b0429e4adfa8cf2e0d02ef9757e8ab96))
* process samples in recalibrate ([a50704c](https://github.com/3d-omics/hg_genotype/commit/a50704c9ab2d41cd1d7b7572075eb728af6f4229))
* process stuff in bed4 file ([42d1e84](https://github.com/3d-omics/hg_genotype/commit/42d1e842d7044889f17f0690b9db2c1cbbec6e11))
* raise slurm latency ([#9](https://github.com/3d-omics/hg_genotype/issues/9)) ([6d1388c](https://github.com/3d-omics/hg_genotype/commit/6d1388c4d1fc48d95254c3b96bad0126e7a346eb))
* recalibrate with cram as input and output ([f718d0e](https://github.com/3d-omics/hg_genotype/commit/f718d0e7401985ec7aa76493b001719844c1a063))
* remove fastp step ([6ae55f0](https://github.com/3d-omics/hg_genotype/commit/6ae55f05b00e8d39a689d52306f8cd293e3e2875))
* Reports by chromosome ([#11](https://github.com/3d-omics/hg_genotype/issues/11)) ([9fe72f1](https://github.com/3d-omics/hg_genotype/commit/9fe72f19e60a7b1bda607af378e7b63bc54b57ec))
* Reports per step and per sample ([#3](https://github.com/3d-omics/hg_genotype/issues/3)) ([a81afe5](https://github.com/3d-omics/hg_genotype/commit/a81afe514316f429a1fae853403361258120e287))
* separate indel and snps, then merge; major rewrite ([5613d8f](https://github.com/3d-omics/hg_genotype/commit/5613d8f5d786c15dd756e033bd8ede6b7436a881))
* set threads to 0 in mark duplicates and apply bsqr since the hard work is done by samtools ([22c838b](https://github.com/3d-omics/hg_genotype/commit/22c838baa5f4efce24625ab4ad259b00581c84a9))
* Set up cache ([#12](https://github.com/3d-omics/hg_genotype/issues/12)) ([e6e2366](https://github.com/3d-omics/hg_genotype/commit/e6e2366a81a81a281c36fd31636e511c3ccf1430))
* sex-dependent ploidy ([#17](https://github.com/3d-omics/hg_genotype/issues/17)) ([295065e](https://github.com/3d-omics/hg_genotype/commit/295065e3be7ec429e00541a53e46bb34918c26ed))
* Slurm-ready ([#5](https://github.com/3d-omics/hg_genotype/issues/5)) ([d7bbd38](https://github.com/3d-omics/hg_genotype/commit/d7bbd38f21332d9f656d93c55fc4f9954d7eabd8))
* substitute picard for gatk4 ([d622d7c](https://github.com/3d-omics/hg_genotype/commit/d622d7c221ecbaf9db0fa72cd88972c344eccf6f))
* update .gitignore ([1c76647](https://github.com/3d-omics/hg_genotype/commit/1c76647b1a73c7613eb3fafa382669895cc45e5c))
* update environments, update pre-commit ([b0de86e](https://github.com/3d-omics/hg_genotype/commit/b0de86e86b33624ce31dcb3b28ad542ea45ae960))
* use bwa-mem2 ([0f29b48](https://github.com/3d-omics/hg_genotype/commit/0f29b48db1e6d54207af66930fa6c53e76a57fb4))


### Bug Fixes

* add environment ([8106036](https://github.com/3d-omics/hg_genotype/commit/8106036a3e05ce4a6d8671c795119bbb58cec8d4))
* add Illumina universal adapters ([baedd26](https://github.com/3d-omics/hg_genotype/commit/baedd266dcf3de21592b2db2f742bca5230d720f))
* add mock bed file ([ab23d03](https://github.com/3d-omics/hg_genotype/commit/ab23d03952930741370873215e96d4ae78f05fb6))
* add mock gtf ([19147c6](https://github.com/3d-omics/hg_genotype/commit/19147c66fa5a7d73dcccc54609bea50763657a19))
* add somalier to env ([3d6f2ed](https://github.com/3d-omics/hg_genotype/commit/3d6f2ed48b57cae02050302d23693de3fbca5599))
* add tbi files ([6a6b8a1](https://github.com/3d-omics/hg_genotype/commit/6a6b8a145ef518918e2f13e274b892e4d765b181))
* add tmp folder for mark_duplicates ([77424d6](https://github.com/3d-omics/hg_genotype/commit/77424d64db153e2b3222502428b3a759b448f265))
* disable level ([cbbd1b5](https://github.com/3d-omics/hg_genotype/commit/cbbd1b5edf1eb2ca3ffd7c9477aa8ac36fbac59b))
* disable max_threads ([aa54952](https://github.com/3d-omics/hg_genotype/commit/aa549523d93b8d17f4a3757560f2047996342a5e))
* disable snpeff ([25b3cbf](https://github.com/3d-omics/hg_genotype/commit/25b3cbf248092614fdd74c2a5b072277a2433c9e))
* do not ignore resources ([ad06ced](https://github.com/3d-omics/hg_genotype/commit/ad06cedfcc5513ad63eaad83a1c34d546cf21e92))
* do not include preprocess ([fc74451](https://github.com/3d-omics/hg_genotype/commit/fc74451f11fc6e58154a367f9b4cca5ba6497a36))
* fastp + gha ([#16](https://github.com/3d-omics/hg_genotype/issues/16)) ([c735460](https://github.com/3d-omics/hg_genotype/commit/c735460ebe2ebd26d70f46bea4fb68f9b987bf98))
* fix downstream side-effects in all rules ([808b49e](https://github.com/3d-omics/hg_genotype/commit/808b49eaf3efa4fb0bf359be662a622317d0e096))
* generate mock vcf file when the chromosome is not present in the sample ([70f6bb5](https://github.com/3d-omics/hg_genotype/commit/70f6bb5e35bb903bd136615be8eeac3d9e1449ee))
* give vep reports a name that multiqc can recognize ([4b9a6d7](https://github.com/3d-omics/hg_genotype/commit/4b9a6d774347bf102be5aa2d3da01b5c50dc7d82))
* I should test more often ([a96251a](https://github.com/3d-omics/hg_genotype/commit/a96251a8a1f8aece22f29c91f0f941b0fe0599e1))
* limit fastp to 16 threads ([#20](https://github.com/3d-omics/hg_genotype/issues/20)) ([8144bcb](https://github.com/3d-omics/hg_genotype/commit/8144bcbee345838f9139165152d12e026cf63a19))
* remove minimizer sort ([506a8de](https://github.com/3d-omics/hg_genotype/commit/506a8de68ea9da0168d43be3324c1674b7e57351))
* safety check to convert fasta to bed4 ([22db278](https://github.com/3d-omics/hg_genotype/commit/22db278057abeac1758b0fd0681090a163466617))
* Slurm times + weekly cache ([#15](https://github.com/3d-omics/hg_genotype/issues/15)) ([5c91a49](https://github.com/3d-omics/hg_genotype/commit/5c91a49e1acd148ea7806e74a93a1e2ee3cb54a3))
* snpeff download db ([c55da24](https://github.com/3d-omics/hg_genotype/commit/c55da24db165a8d194b8eb2d6b22e171f3ca4135))
* sort gtf with bedtools ([fb2eab5](https://github.com/3d-omics/hg_genotype/commit/fb2eab5074b445a60714f5aa4895ad381d6b3c91))
* update documentation ([#19](https://github.com/3d-omics/hg_genotype/issues/19)) ([44a33cd](https://github.com/3d-omics/hg_genotype/commit/44a33cd0c4d3323749a071f6324c2e1bd6430181))


### Performance Improvements

* add somalier find_sites ([dceec42](https://github.com/3d-omics/hg_genotype/commit/dceec4218192d5fc1e81403168a92c6f894359a1))
* use minimal compression for temp file ([111263b](https://github.com/3d-omics/hg_genotype/commit/111263b64f31486662151318389974fb3245c060))
