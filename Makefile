all:
	dmd bamreader.d baminputstream.d bgzfrange.d -ofbamreader -O -release -inline

debug:
	dmd bamreader.d baminputstream.d bgzfrange.d -ofbamreader -debug
