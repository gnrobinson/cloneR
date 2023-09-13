### test-cloneR.R ###
### Note: this unit test tests all functionality of CloneR ###

test_that("Testing CloneR Functionality", {
  subset_num = 5
  snp_num = 10
  Kmin=2
  Kmax=5

  ### Test .vcf files
  ### 1. Ploidy = 2, no sequencing information

  cloneR("data/unit_testing/ploidy2_noinfo.unit_test.vcf", snps = snp_num, subsets = subset_num, K=Kmin:Kmax, plot=FALSE, ploidy = 2)

  # Test input file conversion (VCF)
  concise_file <- read.table("data/unit_testing/ploidy2_noinfo.unit_test.vcf.concise.vcf")
  expect_equal(sum(suppressWarnings(stringr::str_detect(tail(concise_file[,10:ncol(concise_file)]), ':'))), 0)

  # Test genotype file creation
  geno_file <- cloneR::read.geno("data/unit_testing/ploidy2_noinfo.unit_test.vcf.concise.geno")
  expect_s3_class(geno_file, "data.frame")
  expect_equal(ncol(geno_file), 1)
  expect_equal(nrow(geno_file), 8634)

  # Test subsets creation
  if (!dir.exists("subsets")){
    print("subsets directory not successfully created")
  }

  expect_equal(length(list.files(path = "subsets/")), subset_num)

  geno_in <- cloneR::read.geno("subsets/rep_10_1.geno")

  expect_s3_class(geno_in, "data.frame")
  expect_equal(ncol(geno_in), 1)
  expect_equal(nrow(geno_in), snp_num)

  # Test Q files

  expect_equal(length(list.files(path = "q_files/")), length(list.files(path = "subsets/"))*(Kmax-Kmin+1))

  Q_output_K2 <- read.table(file = "q_files/rep_10_1_r1.2.Q")
  Q_output_K3 <- read.table(file = "q_files/rep_10_1_r1.3.Q")
  Q_output_K4 <- read.table(file = "q_files/rep_10_1_r1.4.Q")
  Q_output_K5 <- read.table(file = "q_files/rep_10_1_r1.5.Q")

  expect_equal(ncol(Q_output_K2), 2)
  expect_equal(ncol(Q_output_K3), 3)
  expect_equal(ncol(Q_output_K4), 4)
  expect_equal(ncol(Q_output_K5), 5)

  expect_lt(sum(rowSums(Q_output_K2))/nrow(Q_output_K2), 1.01)
  expect_lt(sum(rowSums(Q_output_K3))/nrow(Q_output_K3), 1.01)
  expect_lt(sum(rowSums(Q_output_K4))/nrow(Q_output_K4), 1.01)
  expect_lt(sum(rowSums(Q_output_K5))/nrow(Q_output_K5), 1.01)

  expect_gt(sum(rowSums(Q_output_K2))/nrow(Q_output_K2), 0.99)
  expect_gt(sum(rowSums(Q_output_K3))/nrow(Q_output_K3), 0.99)
  expect_gt(sum(rowSums(Q_output_K4))/nrow(Q_output_K4), 0.99)
  expect_gt(sum(rowSums(Q_output_K5))/nrow(Q_output_K5), 0.99)

  # Test calc_dist function

  membership_table <- read.csv(file = "membership.txt", sep = '\t')
  membership_rec_table <- read.csv(file = "membership_recommended.txt", sep = '\t')

  expect_s3_class(membership_table, "data.frame")
  expect_equal(ncol(membership_table), Kmax)
  expect_gt(nrow(membership_table), 2)

  expect_s3_class(membership_rec_table, "data.frame")
  expect_equal(ncol(membership_rec_table), 2)
  expect_gt(nrow(membership_rec_table), 2)

  file.remove(c("data/unit_testing/ploidy2_noinfo.unit_test.vcf.concise.geno", "data/unit_testing/ploidy2_noinfo.unit_test.vcf.concise.removed",
                "data/unit_testing/ploidy2_noinfo.unit_test.vcf.concise.vcf", "data/unit_testing/ploidy2_noinfo.unit_test.vcf.concise.vcfsnp",
                "membership.txt", "membership_recommended.txt"))
  unlink("subsets", recursive = TRUE)
  unlink("q_files", recursive = TRUE)

  ### 2. Ploidy = 2, additional sequencing information

  cloneR("data/unit_testing/ploidy2_plusinfo.unit_test.vcf", snps = snp_num, subsets = subset_num, K=Kmin:Kmax, plot=FALSE, ploidy = 2)

  # Test input file conversion (VCF)
  concise_file <- read.table("data/unit_testing/ploidy2_plusinfo.unit_test.vcf.concise.vcf")
  expect_equal(sum(suppressWarnings(stringr::str_detect(tail(concise_file[,10:ncol(concise_file)]), ':'))), 0)

  # Test genotype file creation
  geno_file <- cloneR::read.geno("data/unit_testing/ploidy2_plusinfo.unit_test.vcf.concise.geno")
  expect_s3_class(geno_file, "data.frame")
  expect_equal(ncol(geno_file), 1)
  expect_equal(nrow(geno_file), 9198)

  # Test subsets creation
  if (!dir.exists("subsets")){
    print("subsets directory not successfully created")
  }

  expect_equal(length(list.files(path = "subsets/")), subset_num)

  geno_in <- cloneR::read.geno("subsets/rep_10_1.geno")

  expect_s3_class(geno_in, "data.frame")
  expect_equal(ncol(geno_in), 1)
  expect_equal(nrow(geno_in), snp_num)

  # Test Q files

  expect_equal(length(list.files(path = "q_files/")), length(list.files(path = "subsets/"))*(Kmax-Kmin+1))

  Q_output_K2 <- read.table(file = "q_files/rep_10_1_r1.2.Q")
  Q_output_K3 <- read.table(file = "q_files/rep_10_1_r1.3.Q")
  Q_output_K4 <- read.table(file = "q_files/rep_10_1_r1.4.Q")
  Q_output_K5 <- read.table(file = "q_files/rep_10_1_r1.5.Q")

  expect_equal(ncol(Q_output_K2), 2)
  expect_equal(ncol(Q_output_K3), 3)
  expect_equal(ncol(Q_output_K4), 4)
  expect_equal(ncol(Q_output_K5), 5)

  expect_lt(sum(rowSums(Q_output_K2))/nrow(Q_output_K2), 1.01)
  expect_lt(sum(rowSums(Q_output_K3))/nrow(Q_output_K3), 1.01)
  expect_lt(sum(rowSums(Q_output_K4))/nrow(Q_output_K4), 1.01)
  expect_lt(sum(rowSums(Q_output_K5))/nrow(Q_output_K5), 1.01)

  expect_gt(sum(rowSums(Q_output_K2))/nrow(Q_output_K2), 0.99)
  expect_gt(sum(rowSums(Q_output_K3))/nrow(Q_output_K3), 0.99)
  expect_gt(sum(rowSums(Q_output_K4))/nrow(Q_output_K4), 0.99)
  expect_gt(sum(rowSums(Q_output_K5))/nrow(Q_output_K5), 0.99)

  # Test calc_dist function

  membership_table <- read.csv(file = "membership.txt", sep = '\t')
  membership_rec_table <- read.csv(file = "membership_recommended.txt", sep = '\t')

  expect_s3_class(membership_table, "data.frame")
  expect_equal(ncol(membership_table), Kmax)
  expect_gt(nrow(membership_table), 2)

  expect_s3_class(membership_rec_table, "data.frame")
  expect_equal(ncol(membership_rec_table), 2)
  expect_gt(nrow(membership_rec_table), 2)

  file.remove(c("data/unit_testing/ploidy2_plusinfo.unit_test.vcf.concise.geno", "data/unit_testing/ploidy2_plusinfo.unit_test.vcf.concise.removed",
                "data/unit_testing/ploidy2_plusinfo.unit_test.vcf.concise.vcf", "data/unit_testing/ploidy2_plusinfo.unit_test.vcf.concise.vcfsnp",
                "membership.txt", "membership_recommended.txt"))
  unlink("subsets",recursive = TRUE)
  unlink("q_files",recursive = TRUE)

  ### 3. Ploidy = 1, additional sequencing information

  cloneR("data/unit_testing/ploidy1_plusinfo.unit_test.vcf", snps = snp_num, subsets = subset_num, K=Kmin:Kmax, plot=FALSE, ploidy = 1)

  # Test input file conversion (VCF)
  concise_file <- read.table("data/unit_testing/ploidy1_plusinfo.unit_test.vcf.concise.vcf")
  expect_equal(sum(suppressWarnings(stringr::str_detect(tail(concise_file[,10:ncol(concise_file)]), ':'))), 0)

  # Test genotype file creation
  geno_file <- cloneR::read.geno("data/unit_testing/ploidy1_plusinfo.unit_test.vcf.concise.vcf.geno")
  expect_s3_class(geno_file, "data.frame")
  expect_equal(ncol(geno_file), 1)
  expect_equal(nrow(geno_file), 9912)

  # Test subsets creation
  if (!dir.exists("subsets")){
    print("subsets directory not successfully created")
  }

  expect_equal(length(list.files(path = "subsets/")), subset_num)

  geno_in <- cloneR::read.geno("subsets/rep_10_1.geno")

  expect_s3_class(geno_in, "data.frame")
  expect_equal(ncol(geno_in), 1)
  expect_equal(nrow(geno_in), snp_num)

  # Test Q files

  expect_equal(length(list.files(path = "q_files/")), length(list.files(path = "subsets/"))*(Kmax-Kmin+1))

  Q_output_K2 <- read.table(file = "q_files/rep_10_1_r1.2.Q")
  Q_output_K3 <- read.table(file = "q_files/rep_10_1_r1.3.Q")
  Q_output_K4 <- read.table(file = "q_files/rep_10_1_r1.4.Q")
  Q_output_K5 <- read.table(file = "q_files/rep_10_1_r1.5.Q")

  expect_equal(ncol(Q_output_K2), 2)
  expect_equal(ncol(Q_output_K3), 3)
  expect_equal(ncol(Q_output_K4), 4)
  expect_equal(ncol(Q_output_K5), 5)

  expect_lt(sum(rowSums(Q_output_K2))/nrow(Q_output_K2), 1.01)
  expect_lt(sum(rowSums(Q_output_K3))/nrow(Q_output_K3), 1.01)
  expect_lt(sum(rowSums(Q_output_K4))/nrow(Q_output_K4), 1.01)
  expect_lt(sum(rowSums(Q_output_K5))/nrow(Q_output_K5), 1.01)

  expect_gt(sum(rowSums(Q_output_K2))/nrow(Q_output_K2), 0.99)
  expect_gt(sum(rowSums(Q_output_K3))/nrow(Q_output_K3), 0.99)
  expect_gt(sum(rowSums(Q_output_K4))/nrow(Q_output_K4), 0.99)
  expect_gt(sum(rowSums(Q_output_K5))/nrow(Q_output_K5), 0.99)

  # Test calc_dist function

  membership_table <- read.csv(file = "membership.txt", sep = '\t')
  membership_rec_table <- read.csv(file = "membership_recommended.txt", sep = '\t')

  expect_s3_class(membership_table, "data.frame")
  expect_equal(ncol(membership_table), Kmax)
  expect_gt(nrow(membership_table), 2)

  expect_s3_class(membership_rec_table, "data.frame")
  expect_equal(ncol(membership_rec_table), 2)
  expect_gt(nrow(membership_rec_table), 2)

  file.remove(c("data/unit_testing/ploidy1_plusinfo.unit_test.vcf.concise.vcf.geno",
                "data/unit_testing/ploidy1_plusinfo.unit_test.vcf.concise.vcf",
                "membership.txt", "membership_recommended.txt"))
  unlink("subsets", recursive = TRUE)
  unlink("q_files", recursive = TRUE)
})
