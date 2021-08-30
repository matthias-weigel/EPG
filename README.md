# EPG
The extended phase graph (EPG) concept utilizes a Fourier based magnetization description in terms of so called *configurations states* for depicting and understanding magnetization response in MRI (magnetic resonance imaging) in both a pictorial and quantitative way.
A detailed introduction to the topic including a broader overview of EPGs - considered being helpful by many readers - can be found in [Weigel M. J Magn Reson Imaging 2015. DOI:10.1002/jmri.24619](https://doi.org/10.1002/jmri.24619).


## The EPG codes

The available codes are distributed under the MIT license and essentially belong to the above mentioned manuscript by [M. Weigel](https://doi.org/10.1002/jmri.24619). If you find them useful, please, cite [this repository](https://github.com/matthias-weigel/EPG) and also include a reference to the associated, basic manuscript [Weigel M. J Magn Reson Imaging 2015;41:266-295. DOI:10.1002/jmri.24619](https://doi.org/10.1002/jmri.24619). Thanks a lot!

The example codes represent basic programming examples for multi spin-echo (CP, CPMG, TSE, FSE, RARE, ...)  as well as steady state free precession (SSFP, FLASH, FISP, GRASS, FAST, ...) based MR sequences. Originally, these codes resided at the website [http://epg.matthias-weigel.net](http://epg.matthias-weigel.net) from 2015 to 2018 - until they had to be removed including some additional services due to GDPR issues with the European Union. After that the codes were sent out via email on request.

In the first 3.5 years between 2015 and 2018, the EPG codes were downloaded 715 times from the website, and there is still a continuing interest in the source codes. The author would have *never expected* this success even in his wildest dreams and feels really honored. As a consequence, the original codes are provided in this GitHub repository now.

* **SSFP based sequences: PM-solution** (`ssfp_epg_domain_fplus_fminus.m`)

  Calculates the Extended Phase Graph (EPG) for some variants of gradient echo (GE) / steady state free precession (SSFP) sequences. The particular type of the (idealized) GE/SSFP is specified in the hard coded settings of the software.
  
  For usage, see the code and the accompanying text file. The code uses the Fourier based EPG domains F+(+k), F-(-k), Z(+k), see the [mentioned paper](https://doi.org/10.1002/jmri.24619).
  
  Last modified in 01/2015 (Release Version 2.3).

* **SSFP based sequences: PA-solution** (`ssfp_epg_domain_fplus_alone.m`)

  Calculates the Extended Phase Graph (EPG) for some variants of gradient echo (GE) / steady state free precession (SSFP) sequences. The particular type of the (idealized) GE/SSFP is specified in the hard coded settings of the software.
  
  For usage, see the code and the accompanying text file. The code uses the Fourier based EPG domains F+(+k), F+(-k), Z(+k), see the [mentioned paper](https://doi.org/10.1002/jmri.24619).
  
  Last modified in 01/2015 (Release Version 1.8).


* **Multi spin-echo sequences of type CP or CPMG** (`cp_cpmg_epg_domain_fplus_fminus.m`)

  Calculates the Extended Phase Graph (EPG) over one TR / one echo train for multi spin-echo (multi SE) sequences obeying the Carr-Purcell (CP) or Carr-Purcell-Meiboom-Gill (CPMG) conditions. The particular type (CP/CPMG) of the (idealized) spin echo sequence is specified in the hard coded settings of the software. Common CPMG based imaging sequences have acronyms such as SE, TSE, FSE or RARE.
  
  For usage, see the code and the accompanying text file. The code uses the Fourier based EPG domains F+(+k), F-(-k), Z(+k), see the [mentioned paper](https://doi.org/10.1002/jmri.24619).
  
  Last modified in 07/2015 (Release Version 1.3).

## ***!!! HAVE FUN !!!***
