/*
*  Test input fast5
*/

process testInputFast5 {
    tag { fast5 }
        
    input:
    path(fast5) 

    output:
    stdout()

    script:
    """
    fast5_type.py ${fast5}
    """
}

 
 workflow TEST_FAST5 {
    take: fast5
    main:
        testInputFast5(fast5.first())
    	fast5_type = testInputFast5.out.map {  it.trim().toInteger() }
	  	fast5_type.map{it == 0 ? "Single Fast5 files detected!": "MultiFast5 files detected!" }.view()

   emit:
        fast5_type
}

