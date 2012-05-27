require './bindings/libbam.rb'
require 'rspec/expectations'

Before do
  @bam = BamFile.new 'test/data/ex1_header.bam'
end

When /^I use its 'alignments' method$/ do
  @bam.should respond_to(:alignments)
end

Then /^I should be able to iterate the returned object with 'each'$/ do
  @bam.alignments.should respond_to(:each)
end

Then /^the objects which I iterate over should represent the alignments$/ do
  @bam.alignments.each do |read|
    read.should be_instance_of(Alignment)
  end
end

Then /^I should be able to access all fields mentioned in SAM\/BAM format specification$/ do
  @bam.rewind!
  @read = @bam.alignments.first
  @read.read_name.should == 'EAS56_57:6:190:289:82'
  @read.sequence.should == 'CTCAAGGTTGTTGCAAGGGGGTCTATGTGAACAAA'
  @read.position.should == 99
  @read.flag.should == 69
  @read.mapping_quality.should == 0
  @read.cigar_string.should == ''
  @bam.reference_sequences[@read.ref_id].name.should == 'chr1'
  @read.quality.should == [27, 27, 27, 22, 27, 27, 27, 26, 27, 27, 27, 27, 27, 27, 27, 27, 23, 26, 26, 27, 22, 26, 19, 27, 26, 27, 26, 26, 26, 26, 26, 24, 19, 27, 26]
end

Given /^I have an alignment$/ do
  @read = @bam.alignments.first
end

Given /^it contains some tags$/ do
end

When /^I access it like a hash$/ do
  @read.should respond_to(:[])
end

When /^I use 2-character string as a key$/ do
  @key = 'MF'
end

When /^the alignment has such tag$/ do
  @read[@key].should_not be_nil
end

Then /^I should be able to see corresponding value$/ do
  @read[@key].should be == 192
end

Then /^it should be a simple Ruby object \(Array, Numeric, or String\)$/ do
  @read[@key].should be_kind_of Numeric
end

When /^I use string of length different than two, as a key,$/ do
  @key = 'key'
end

Then /^exception should be thrown\.$/ do
  expect{@read[@key]}.to raise_error(RuntimeError)
end

When /^it doesn't contain the requested tag$/ do
  @key = 'hq'
end

Then /^nil should be returned\.$/ do
  @read[@key].should be_nil
end

When /^I use its 'tags' method$/ do
  @tags = @read.tags
end

Then /^I should be able to work with the returned object just like with Hash$/ do
  @tags.should be_kind_of Hash
  @tags['MF'].should be == 192
  @tags.keys.should be == ['MF']
  @tags.values.should be == [192]
end
