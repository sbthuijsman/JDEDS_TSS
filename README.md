# JDEDS_TSS
These are all the files used to compute the results for "Transformational Supervisor Synthesis for Evolving Systems", Thuijsman S., Reniers M. (2022), Journal of Discrete Event Dynamic Systems

The implementation has been corrected after publication. This is the corrected implementation. A correction notice has been submitted to the journal.

The experiments can be repeated by performing:
- MLsynth.m (for base and variant model); this is the supervisor synthesis algorithm (Algorithm 1)
- MLdeltasynth.m (for base model, variant model, and modeldelta); performs transformational supervisor synthesis. doITSS=true will perform ITSS, doITSS=false will perform GTSS.

The models are found in the MLsys folder. 
The reachable variant of a model is indicated by '_r'.
