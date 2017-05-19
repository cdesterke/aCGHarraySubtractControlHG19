# aCGHarraySubtractControlHG19
Shell script which allowed to remove european polymorphism in CNV and annotate genomic intervals on HG19 after subtracting genomic intervals of control condition

This script needs the installation of Bedtools utilities in the OS environement. This script also needs dependencies: installing 5 annotation databases in the subfolder DATABASES_HG19:

   - hg19GENEname_sorted.bed: Refseq HG19 annotations
   - scandbHG19CEU_sorted.bed: CNV polymorphisms in european population, data derived from http://www.scandb.org/ website
   - plurinetHG19_sorted.bed: genes that are characterized as belonging to a pluripotency network (PLURINET DB)
   - CosmicHG19_sorted.bed: genes presenting somatic alteration in cancer and leukemia
   - TFHG19.bed: census database of molecules classed as transcription factor

Script needs to be execute in BASH environnement with previous dependencies.

USAGE: 
      
      $sh acgh_subcontrol_HG19.sh example.bed control.bed

EXAMPLE FILES: "example.bed" and "control.bed" with minimum 5 columns needs to be pass as first parameter of the shell line bed file and control bed file as 2nd parameter of shell line, bed files need to be place at the parental directory so same the directory than the script! 

DATABASES :

   - Transcription factor database (Vaquerizas et al., 2009)
   - COSMIC database (Futreal et al., 2004)
   - Plurinet database (Müller et al., 2008)
   - European CNV polymorphism annotation (Gamazon et al., 2010)
   - RefSeq annotation database on HG18 (O’Leary et al., 2016)

DEPENDENCIES:

   - Bedtools (Quinlan, 2014)

References

-Futreal, P.A., Coin, L., Marshall, M., Down, T., Hubbard, T., Wooster, R., Rahman, N., and Stratton, M.R. (2004). A census of human cancer genes. Nat. Rev. Cancer 4, 177–183.

-Gamazon, E.R., Zhang, W., Konkashbaev, A., Duan, S., Kistner, E.O., Nicolae, D.L., Dolan, M.E., and Cox, N.J. (2010). SCAN: SNP and copy number annotation. Bioinforma. Oxf. Engl. 26, 259–262.

-Müller, F.-J., Laurent, L.C., Kostka, D., Ulitsky, I., Williams, R., Lu, C., Park, I.-H., Rao, M.S., Shamir, R., Schwartz, P.H., et al. (2008). Regulatory networks define phenotypic classes of human stem cell lines. Nature 455, 401–405.

-O’Leary, N.A., Wright, M.W., Brister, J.R., Ciufo, S., Haddad, D., McVeigh, R., Rajput, B., Robbertse, B., Smith-White, B., Ako-Adjei, D., et al. (2016). Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. Nucleic Acids Res. 44, D733-745.

-Quinlan, A.R. (2014). BEDTools: The Swiss-Army Tool for Genome Feature Analysis. Curr. Protoc. Bioinforma. 47, 11.12.1-34.

-Vaquerizas, J.M., Kummerfeld, S.K., Teichmann, S.A., and Luscombe, N.M. (2009). A census of human transcription factors: function, expression and evolution. Nat. Rev. Genet. 10, 252–263.
