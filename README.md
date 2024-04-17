# inclXSEC

To use it:
1) unzip it
2) run "make"
3) ./distrib_x86_64 < "name of input file" or just ./distrib_x86_64 but you will need to manually type all the parameters. Finally, you can create a script similar to my runOsp.py to run the code for multiple input parameters.
4) the output file and kinematics are specified in input parameters.

A few comments:
a) distrib.f is the file were input and output specified (line 47): write(33,20) w,(theta*raddeg),p_el,sig,pureel,elra,eltail,(sig_rc-elra-eltail),sig_rc,pureel_noWidth

so the output file is:
W, theta in deg, momentum GeV,sig (inelastic non RAD XSEC), pureel (ignore I need it for my studies), elra (elastic XSEC with radiation ON), eltail (elastic tail XSEC), sig_rc-elra-eltail (inelastic RAD XSEC), sig_rc (full RAD XSEC eles rad + tail + inelastic rad), pureel_noWidth (ignore it)

so you can modify line 47 if you want a different output.

b) After any change in the code you need to run "make" again. You do not need 'make' when you change input parameters only.

c) External RC after scattering (straggling after the interaction point) are OFF because GEANT takes care of it.
