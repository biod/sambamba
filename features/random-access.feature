Feature: random access to BAM file
    In order to retrieve information about specific regions,
    I want to be able to quickly fetch alignments overlapping a region.

    Scenario: fetching alignments
        Given I have a BAM file
          And it's sorted by coordinate
          And I have its index as well
         When I specify reference sequence and region (0-based beginning and end positions)
         Then I should be able to immediately have access to alignments overlapping it
