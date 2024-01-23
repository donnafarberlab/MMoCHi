# Landmark Registration Advice

Manually checking landmark registration can be time consuming, but it's a useful way to ensure effective ADT integration. One tricky step of landmark registration is setting effective `marker_bandwidth` values to maximize on automatic peak detection. In our experience with multiple CITE-seq panels (TotalSeq A, B, and C), the best marker_bandwidth value for individual antibodies are consistent across panels. Thus, to help jumpstart the process, below is a list of some of our modifications to `marker_bandwidth` for various antibodies. For markers not listed, the built-in default (`0.2`) seems to work consistently.

We hope to continue to update this page, expanding advice to additional antibody panels and sample types! If you have any input, click [here](https://github.com/donnafarberlab/MMoCHi/issues/), and create a "New issue" button, choose "Submit advice".

## Bandwidths for human immune cell markers 
(Applies to the Universal TotalSeq-A, B, and C panels)

| Marker Name | marker_bandwidth  |
|         ---:|      :---:        |
| CD2 | 0.1 |
| CD3 | 0.1 |
| CD7 | 0.1 |
| CD16 | 0.15 | 
| CD19 | 0.15 | 
| CD28 | 0.5 | 
| CD32 | 0.15 | 
| CD35 | 0.175 | 
| CD40 | 0.5 | 
| CD56 (NCAM) | 0.15 | 
| CD62L | 0.1 | 
| CD82 | 0.3 | 
| CD94 | 0.15 | 
| CD107a (LAMP-1) | 0.5 | 
| CD123 | 0.6 | 
| CD127 (IL-7Ra) | 0.15 |
| CD163 | 0.15 |
| CD185 (CXCR5) | 0.3 | 
| CD335 (NKp46) | 0.15 | 
| Ig light chain k | 0.3 | 

