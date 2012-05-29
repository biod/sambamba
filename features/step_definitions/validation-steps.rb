Before do
  # this file is known to contain some invalid records
  @tagsbam = BamFile.new './test/data/tags.bam'
end

Given /^I have an alignment from a BAM file$/ do
  @alignment = @tagsbam.alignments.to_a[32]
end

When /^I call 'valid\?' method$/ do
  @is_valid = @alignment.valid?
end

Then /^it should return whether it is valid or not$/ do
  @is_valid.should be_true
end

Given /^I have a BAM file$/ do
  @tagsbam.rewind!
end

When /^I want to iterate over its records$/ do
  @records = @tagsbam.alignments
end

Then /^I should have an option to skip invalid ones$/ do
  @records.should respond_to(:each).with(1).argument
end

Then /^all the reads in this case should be valid$/ do
  @records.each(:valid => true) do |record|
    record.should be_valid
  end
end
