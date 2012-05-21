Feature: iterating alignment records

    In order to have access to all information contained in a BAM file,
    As a bioinformatician,
    I want to be able to iterate alignment records from Ruby
    And have access to all their fields and tags.

    Scenario: accessing alignment records
        Given I opened a valid BAM file
        When I use its 'alignments' method
        Then I should be able to iterate the returned object with 'each'
        And the objects which I iterate over should represent the alignments
        And I should be able to access all fields mentioned in SAM/BAM format specification

    Scenario: access existing alignment tag
        Given I have an alignment
          And it contains some tags
         When I access it like a hash
          And I use 2-character string as a key
          And the alignment has such tag
         Then I should be able to see corresponding value
          And it should be a simple Ruby object (Array, Numeric, or String)

    Scenario: invalid tag key (not of length 2)
        Given I have an alignment
         When I access it like a hash
          But I use string of length different than two, as a key,
         Then exception should be thrown.

    Scenario: accessing non-existing alignment tag
        Given I have an alignment
          And it contains some tags
         When I access it like a hash
          But it doesn't contain the requested tag
          Then nil should be returned.

    Scenario: fetching all tags as a hash
        Given I have an alignment
         When I use its 'tags' method
         Then I should be able to work with the returned object just like with Hash
