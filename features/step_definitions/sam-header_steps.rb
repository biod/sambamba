Given /^I opened a valid BAM file$/ do
  filename = './test/data/ex1_header.bam'
  File.exists?(filename).should be_true
  @bamfile = BamFile.new filename
end

Given /^it contains SAM header$/ do
  @bamfile.header.raw_contents.length.should be > 0
end

When /^I call 'header' method$/ do
  @header = @bamfile.header
end

Then /^I should see text of SAM header$/ do
  @header.raw_contents.should be_kind_of String
end

Given /^SAM header contains @HD line$/ do
  @header = @bamfile.header
  @header.raw_contents.should =~ /^@HD/
end

Then /^I should be able to see format version$/ do
  @version = @header.version
  @version.should be_kind_of String
  @version.length.should be > 0
end

Then /^I should be able to see sorting order$/ do
  @sorting_order = @header.sorting_order
  @sorting_order.should be_kind_of String
  @sorting_order.length.should be > 0
end

Given /^SAM header contains @SQ lines$/ do
  @header = @bamfile.header
  @header.sq_lines.length.should be > 0
end

Then /^I should be able to iterate them$/ do
  @sq_lines = @header.sq_lines
  @sq_lines.should be_kind_of Array
end

Then /^I should be able to see sequence names$/ do
  @line = @sq_lines.first
  @line.should respond_to(:sequence_name).with(0).arguments
  @line.sequence_name.should be_kind_of String
  @line.sequence_name.length.should be > 0
end

Then /^I should be able to see their lengths$/ do
  @line.should respond_to(:sequence_length).with(0).arguments
  @line.sequence_length.should be_kind_of Numeric
end
