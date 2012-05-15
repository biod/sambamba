Feature: access to information from SAM header

    In order to work with BAM file,
    I want to see what its header contains.

    Background:
        Given I opened a valid BAM file
          And it contains SAM header

    Scenario: getting raw text
         When I call 'header' method
         Then I should see text of SAM header

    Scenario: accessing version and sorting order
         When SAM header contains @HD line
         Then I should be able to see format version
          And I should be able to see sorting order

    Scenario: getting information about reference sequences
         When SAM header contains @SQ lines
         Then I should be able to iterate them
          And I should be able to see sequence names
          And I should be able to see their lengths
