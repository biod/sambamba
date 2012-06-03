Before do
    @bam = BamFile.new './test/data/ex1_header.bam'
end

Given /^it's sorted by coordinate$/ do
    @bam.header.sorting_order.should == 'coordinate'
end

Given /^I have its index as well$/ do
    @bam.should have_index
end

When /^I specify reference sequence and region \(0-based beginning and end positions\)$/ do
    @region = (1400 ... 1500)
    @chr = "chr2"
end

Then /^I should be able to immediately have access to alignments overlapping it$/ do
    @alignments = @bam.fetch @chr, @region
    @alignments.should respond_to(:each).with(0).arguments
    @alignments.to_a.length.should == 75
end
