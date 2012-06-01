Feature: syntax sugar

    In order to enjoy writing my scripts,
    As a Rubyista,
    I want some syntax sugar.

    Scenario: fetching alignments
        Given I have a BAM file
          And associated BAI file
         When I say "bam.alignments.referencing(chromosome).overlapping(500.kbp .. 600.kbp)"
         Then I should get these alignments

    Scenario: using shortcuts
        Given I have a BAM file
          And associated BAI file
         When I say "bam[chromosome][500.kbp .. 600.kbp]"
         Then I should get these alignments
