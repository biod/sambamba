Feature: alignment validation

    In order to be able to filter out invalid reads,
    As a developer,
    I want validation support.

    Scenario: checking single read
        Given I have an alignment from a BAM file
         When I call 'valid?' method
         Then it should return whether it is valid or not

    Scenario: iterating over valid records
        Given I have a BAM file
         When I want to iterate over its records
         Then I should have an option to skip invalid ones
          And all the reads in this case should be valid
