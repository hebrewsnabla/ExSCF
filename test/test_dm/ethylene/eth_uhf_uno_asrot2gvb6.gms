----- GAMESS execution script 'rungms' -----
This job is running on host xn01
under operating system Linux at Wed Jan 20 01:39:19 CST 2021
Available scratch disk space (Kbyte units) at beginning of the job is
Filesystem     1K-blocks  Used Available Use% Mounted on
/dev/sda5      150914052 54728 150859324   1% /scratch
GAMESS temporary binary files will be written to /scratch/scr/srwang/gamess
GAMESS supplementary output files will be written to /scratch/scr/srwang/gamess
Copying input file ./eth_uhf_uno_asrot2gvb6.inp to your run's scratch directory...
cp ./eth_uhf_uno_asrot2gvb6.inp /scratch/scr/srwang/gamess/eth_uhf_uno_asrot2gvb6.F05
unset echo
/share/home/srwang/gamess/ddikick.x /share/home/srwang/gamess/gamess.01.x eth_uhf_uno_asrot2gvb6 -ddi 1 1 xn01 -scr /scratch/scr/srwang/gamess

 Distributed Data Interface kickoff program.
 Initiating 1 compute processes on 1 nodes to run the following command:
 /share/home/srwang/gamess/gamess.01.x eth_uhf_uno_asrot2gvb6 

          *******************************************************
          *                                                     *
          *          GAMESS VERSION = 30 JUN 2020 (R1)          *
          *                                                     *
          *              doi.org/10.1063/5.0005188              *
          *                                                     *
          **************** 64 BIT INTEL VERSION *****************

  GAMESS HAS BEEN MADE POSSIBLE WITH IMPORTANT CONTRIBUTIONS FROM
  THE FOLLOWING INDIVIDUALS (IN ALPHABETICAL ORDER):

     IVANA ADAMOVIC, CHRISTINE AIKENS, TOMOKO AKAMA,
     YURI ALEXEEV, YURIKO AOKI, POOJA ARORA, TOSHIO ASADA,
     ANDREY ASADCHEV, KIM K. BALDRIDGE, PRADIPTA BANDYOPADHYAY,
     MARIA BARYSZ, ROB BELL, JONATHAN BENTZ, COLLEEN BERTONI,
     JERRY A. BOATZ, BRETT BODE, KURT BRORSEN, KURT R. BRORSEN,
     LAIMUTIS BYTAUTAS, CALEB CARLIN, LAURA CARRINGTON,
     GALINA CHABAN, BENOIT CHAMPAGNE, WEI CHEN, MAHITO CHIBA,
     DAN CHIPMAN, CHEOL HO CHOI, TANNER CULPITT, PAUL DAY,
     ALBERT DEFUSCO, NUWAN DESILVA, J. EMILIANO DEUSTUA,
     TIM DUDLEY, MICHEL DUPUIS, STEVEN T. ELBERT,
     DMITRI FEDOROV, ALEX FINDLATER, GRAHAM FLETCHER,
     MARK FREITAG, CHRISTIAN FRIEDL, DAVID GARMER,
     IGOR S. GERASIMOV, KURT GLAESEMANN, MARK S. GORDON,
     JEFFREY GOUR, FENG LONG GU, EMILIE B. GUIDEZ,
     ANASTASIA GUNINA, SHARON HAMMES-SCHIFFER, MASATAKE HARUTA,
     KIMIHIKO HIRAO, YASUHIRO IKABATA, TZVETELIN IORDANOV,
     STEPHAN IRLE, KAZUYA ISHIMURA, JOE IVANIC, FRANK JENSEN,
     JAN H. JENSEN, VISVALDAS KAIRYS, MUNEAKI KAMIYA,
     MICHIO KATOUDA, NAOAKI KAWAKAMI, DAN KEMP, BERNARD KIRTMAN,
     KAZUO KITAURA, MARIUSZ KLOBUKOWSKI, MASATO KOBAYASHI,
     PRAKASHAN KORAMBATH, JACEK KORCHOWIEC, SHIRO KOSEKI,
     KAROL KOWALSKI, JIMMY KROMANN, STANISLAW KUCHARSKI,
     HENRY KURTZ, SAROM SOK LEANG, HUI LI, SHUHUA LI, WEI LI,
     JESSE LUTZ, ALEKSANDR O. LYKHIN, MARCIN MAKOWSKI,
     JOANI MATO, NIKITA MATSUNAGA, BENEDETTA MENNUCCI,
     GRANT MERRILL, NORIYUKI MINEZAWA, VLADIMIR MIRONOV,
     EISAKU MIYOSHI, JOHN A. MONTGOMERY JR., HIROTOSHI MORI,
     JONATHAN MULLIN, MONIKA MUSIAL, SHIGERU NAGASE,
     TAKESHI NAGATA, HIROMI NAKAI, TAKAHITO NAKAJIMA,
     YUYA NAKAJIMA, HARUYUKI NAKANO, HIROYA NAKATA, SEAN NEDD,
     HEATHER NETZLOFF, KIET A. NGUYEN, YOSHIO NISHIMOTO,
     BOSILJKA NJEGIC, RYAN OLSON, MIKE PAK, ROBERTO PEVERATI,
     BUU PHAM, PIOTR PIECUCH, ANNA POMOGAEVA, DAVID POOLE,
     SPENCER PRUITT, OLIVIER QUINET, LUKE ROSKOP,
     KLAUS RUEDENBERG, ANDREW SAND, TOSAPORN SATTASATHUCHANA,
     NOZOMI SAWADA, MICHAEL W. SCHMIDT, PATRICK E. SCHNEIDER,
     JUNJI SEINO, PRACHI SHARMA, JUN SHEN, JIM SHOEMAKER,
     YINAN SHU, DEJUN SI, JONATHAN SKONE, LYUDMILA SLIPCHENKO,
     TONY SMITH, JIE SONG, MARK SPACKMAN, CASPER STEINMANN,
     WALT STEVENS, PEIFENG SU, SHUJUN SU, CHET SWALINA,
     TETSUYA TAKETSUGU, ZHEN TAO, NANDUN THELLAMUREGE
     SEIKEN TOKURA, JACOPO TOMASI, TSUGUKI TOUMA, TAKAO TSUNEDA,
     HIROAKI UMEDA, JORGE LUIS GALVEZ VALLEJO, YALI WANG,
     SIMON WEBB, AARON WEST, THERESA L. WINDUS, MARTA WLOCH,
     PENG XU, KIYOSHI YAGI, SUSUMU YANAGISAWA, YANG YANG,
     SOOHAENG YOO, TAKESHI YOSHIKAWA, FEDERICO ZAHARIEV,
     TOBY ZENG

  WHO ARE SUPPORTED BY THEIR INSTITUTION/UNIVERSITY/COMPANY/GROUP
  (IN ALPHABETICAL ORDER):

     EP ANALYTICS, FACULTES UNIVERSITAIRES NOTRE-DAME DE LA PAIX,
     INSTITUTE FOR MOLECULAR SCIENCE, IOWA STATE UNIVERSITY,
     JACKSON STATE UNIVERSITY, JOHANNES KEPLER UNIVERSITY LINZ,
     KYUSHU UNIVERSITY, MICHIGAN STATE UNIVERSITY,
     MIE UNIVERSITY, MOSCOW STATE UNIVERSITY,
     N. COPERNICUS UNIVERSITY, NAGOYA UNIVERSITY,
     NANJING UNIVERSITY, NAT. INST. OF ADVANCED INDUSTRIAL
     SCIENCE AND TECHNOLOGY, NATIONAL INST. OF STANDARDS AND
     TECHNOLOGY, NESMEYANOV INSTITUTE OF ORGANOELEMENT COMPOUNDS
     OF RUSSIAN ACADEMY OF SCIENCES, OSAKA PREFECTURE UNIVERSITY,
     PENNSYLVANIA STATE UNIVERSITY, TOKYO INSTITUTE OF TECHNOLOGY,
     UNIVERSITY OF AARHUS, UNIVERSITY OF ALBERTA,
     UNIVERSITY OF CALIFORNIA AT SANTA BARBARA,
     UNIVERSITY OF COPENHAGEN, UNIVERSITY OF IOWA,
     UNIVERSITY OF MEMPHIS, UNIVERSITY OF MINNESOTA,
     UNIVERSITY OF NEBRASKA, UNIVERSITY OF NEW ENGLAND,
     UNIVERSITY OF NOTRE DAME, UNIVERSITY OF PISA,
     UNIVERSITY OF SILESIA, UNIVERSITY OF TOKYO,
     UNIVERSITY OF ZURICH, WASEDA UNIVERSITY, YALE UNIVERSITY

  GAMESS SOFTWARE MANAGEMENT TEAM FOR THIS RELEASE:

     MARK S. GORDON (TEAM LEAD),
     BRETT BODE (SENIOR ADVISOR),
     GIUSEPPE BARCA,
     COLLEEN BERTONI,
     KRISTOPHER KEIPERT (ADVISOR),
     SAROM S. LEANG (DEVELOPMENT LEAD),
     BUU PHAM,
     JORGE LUIS GALVEZ VALLEJO (WEBSITE ADMINISTRATOR),
     PENG XU

 EXECUTION OF GAMESS BEGUN Wed Jan 20 01:39:19 2021

            ECHO OF THE FIRST FEW INPUT CARDS -
 INPUT CARD> $CONTRL SCFTYP=GVB RUNTYP=ENERGY ICHARG=0 MULT=1 NOSYM=1 ICUT=10               
 INPUT CARD>  MAXIT=500 RELWFN=DK ISPHER=1 $END                                             
 INPUT CARD> $SYSTEM MWORDS=125 $END                                                        
 INPUT CARD> $RELWFN NORDER=1 $END                                                          
 INPUT CARD> $SCF NCO=2 NPAIR=6 DIRSCF=.TRUE. $END                                          
 INPUT CARD> $GUESS GUESS=MOREAD NORB=48 $END                                               
 INPUT CARD> $DATA                                                                          
 INPUT CARD>GAMESS inp format file produced by MOKIT, nbf=48                                
 INPUT CARD>C1   1                                                                          
 INPUT CARD>C     6.         0.73236300        -0.00001200        -0.00002000               
 INPUT CARD>   S  5                                                                         
 INPUT CARD>   1   1.23840169E+03   5.51736505E-03   0.00000000E+00                         
 INPUT CARD>   2   1.86290050E+02   4.10888286E-02   0.00000000E+00                         
 INPUT CARD>   3   4.22511763E+01   1.82253821E-01   0.00000000E+00                         
 INPUT CARD>   4   1.16765579E+01   4.68284594E-01   0.00000000E+00                         
 INPUT CARD>   5   3.59305065E+00   4.45758173E-01   0.00000000E+00                         
 INPUT CARD>   S  1                                                                         
 INPUT CARD>   1   4.02451474E-01   1.00000000E+00   0.00000000E+00                         
 INPUT CARD>   S  1                                                                         
 INPUT CARD>   1   1.30901827E-01   1.00000000E+00   0.00000000E+00                         
 INPUT CARD>   P  3                                                                         
 INPUT CARD>   1   9.46809706E+00   5.68883320E-02   0.00000000E+00                         
 INPUT CARD>   2   2.01035451E+00   3.12940593E-01   0.00000000E+00                         
 INPUT CARD>   3   5.47710047E-01   7.60650165E-01   0.00000000E+00                         
 INPUT CARD>   P  1                                                                         
 INPUT CARD>   1   1.52686138E-01   1.00000000E+00   0.00000000E+00                         
 INPUT CARD>   D  1                                                                         
 INPUT CARD>   1   8.00000000E-01   1.00000000E+00   0.00000000E+00                         
 INPUT CARD>                                                                                
 INPUT CARD>H     1.         1.30178800        -0.65477800         0.65483500               
 INPUT CARD>   S  3                                                                         
 INPUT CARD>   1   1.30107010E+01   3.34854848E-02   0.00000000E+00                         
 INPUT CARD>   2   1.96225720E+00   2.34721871E-01   0.00000000E+00                         
 INPUT CARD>   3   4.44537960E-01   8.13770285E-01   0.00000000E+00                         
 INPUT CARD>   S  1                                                                         
 INPUT CARD>   1   1.21949620E-01   1.00000000E+00   0.00000000E+00                         
 INPUT CARD>   P  1                                                                         
 INPUT CARD>   1   8.00000000E-01   1.00000000E+00   0.00000000E+00                         
 INPUT CARD>                                                                                
 INPUT CARD>H     1.         1.30180201         0.65484800        -0.65476200               
 INPUT CARD>   S  3                                                                         
 INPUT CARD>   1   1.30107010E+01   3.34854848E-02   0.00000000E+00                         
 INPUT CARD>   2   1.96225720E+00   2.34721871E-01   0.00000000E+00                         
 INPUT CARD>   3   4.44537960E-01   8.13770285E-01   0.00000000E+00                         
 INPUT CARD>   S  1                                                                         
 INPUT CARD>   1   1.21949620E-01   1.00000000E+00   0.00000000E+00                         
 INPUT CARD>   P  1                                                                         
 INPUT CARD>   1   8.00000000E-01   1.00000000E+00   0.00000000E+00                         
 INPUT CARD>                                                                                
 INPUT CARD>C     6.        -0.73236300        -0.00001200         0.00002000               
  125000000 WORDS OF MEMORY AVAILABLE


     RUN TITLE
     ---------
 GAMESS inp format file produced by MOKIT, nbf=48                                

 THE POINT GROUP OF THE MOLECULE IS C1      
 THE ORDER OF THE PRINCIPAL AXIS IS     1

 ATOM      ATOMIC                      COORDINATES (BOHR)
           CHARGE         X                   Y                   Z
 C           6.0     1.3839653935       -0.0000226767       -0.0000377945
 H           1.0     2.4600226141       -1.2373510028        1.2374587172
 H           1.0     2.4600490892        1.2374832836       -1.2373207672
 C           6.0    -1.3839653935       -0.0000226767        0.0000377945
 H           1.0    -2.4600226141       -1.2373510028       -1.2374606069
 H           1.0    -2.4600490892        1.2374832836        1.2373207672

          INTERNUCLEAR DISTANCES (ANGS.)
          ------------------------------

                1 C          2 H          3 H          4 C          5 H     

   1 C       0.0000000    1.0871055 *  1.0871014 *  1.4647260 *  2.2350107 *
   2 H       1.0871055 *  0.0000000    1.8520703 *  2.2350104 *  2.9144204 *
   3 H       1.0871014 *  1.8520703 *  0.0000000    2.2350410 *  2.9144127 *
   4 C       1.4647260 *  2.2350104 *  2.2350410 *  0.0000000    1.0871061 *
   5 H       2.2350107 *  2.9144204 *  2.9144127 *  1.0871061 *  0.0000000  
   6 H       2.2350410 *  2.9144127 *  2.9143794 *  1.0871014 *  1.8520711 *

                6 H     

   1 C       2.2350410 *
   2 H       2.9144127 *
   3 H       2.9143794 *
   4 C       1.0871014 *
   5 H       1.8520711 *
   6 H       0.0000000  

  * ... LESS THAN  3.000


     ATOMIC BASIS SET
     ----------------
 THE CONTRACTED PRIMITIVE FUNCTIONS HAVE BEEN UNNORMALIZED
 THE CONTRACTED BASIS FUNCTIONS ARE NOW NORMALIZED TO UNITY

  SHELL TYPE  PRIMITIVE        EXPONENT          CONTRACTION COEFFICIENT(S)

 C         

      1   S       1          1238.4016900    0.005517365054
      1   S       2           186.2900500    0.041088828628
      1   S       3            42.2511763    0.182253821124
      1   S       4            11.6765579    0.468284594318
      1   S       5             3.5930506    0.445758173303

      2   S       6             0.4024515    1.000000000000

      3   S       7             0.1309018    1.000000000000

      4   P       8             9.4680971    0.056888332011
      4   P       9             2.0103545    0.312940593061
      4   P      10             0.5477100    0.760650165148

      5   P      11             0.1526861    1.000000000000

      6   D      12             0.8000000    1.000000000000

 H         

      7   S      13            13.0107010    0.033485484805
      7   S      14             1.9622572    0.234721871038
      7   S      15             0.4445380    0.813770285133

      8   S      16             0.1219496    1.000000000000

      9   P      17             0.8000000    1.000000000000

 H         

     10   S      18            13.0107010    0.033485484805
     10   S      19             1.9622572    0.234721871038
     10   S      20             0.4445380    0.813770285133

     11   S      21             0.1219496    1.000000000000

     12   P      22             0.8000000    1.000000000000

 C         

     13   S      23          1238.4016900    0.005517365054
     13   S      24           186.2900500    0.041088828628
     13   S      25            42.2511763    0.182253821124
     13   S      26            11.6765579    0.468284594318
     13   S      27             3.5930506    0.445758173303

     14   S      28             0.4024515    1.000000000000

     15   S      29             0.1309018    1.000000000000

     16   P      30             9.4680971    0.056888332011
     16   P      31             2.0103545    0.312940593061
     16   P      32             0.5477100    0.760650165148

     17   P      33             0.1526861    1.000000000000

     18   D      34             0.8000000    1.000000000000

 H         

     19   S      35            13.0107010    0.033485484805
     19   S      36             1.9622572    0.234721871038
     19   S      37             0.4445380    0.813770285133

     20   S      38             0.1219496    1.000000000000

     21   P      39             0.8000000    1.000000000000

 H         

     22   S      40            13.0107010    0.033485484805
     22   S      41             1.9622572    0.234721871038
     22   S      42             0.4445380    0.813770285133

     23   S      43             0.1219496    1.000000000000

     24   P      44             0.8000000    1.000000000000

 TOTAL NUMBER OF BASIS SET SHELLS             =   24
 NUMBER OF CARTESIAN GAUSSIAN BASIS FUNCTIONS =   50
 NOTE: THIS RUN WILL RESTRICT THE MO VARIATION SPACE TO SPHERICAL HARMONICS.
 THE NUMBER OF ORBITALS KEPT IN THE VARIATIONAL SPACE WILL BE PRINTED LATER.
 NUMBER OF ELECTRONS                          =   16
 CHARGE OF MOLECULE                           =    0
 SPIN MULTIPLICITY                            =    1
 NUMBER OF OCCUPIED ORBITALS (ALPHA)          =    8
 NUMBER OF OCCUPIED ORBITALS (BETA )          =    8
 TOTAL NUMBER OF ATOMS                        =    6
 THE NUCLEAR REPULSION ENERGY IS       31.6688684670

     ---------------------------------
     RELATIVISTIC WAVEFUNCTION OPTIONS
     ---------------------------------
     METHOD=DK       NORDER= 1       NESOC=  0      MODEQR=      1
     NRATOM=  0      QMTTOL= 1.000E-10      QRTOL= 1.000E-02
     CLIGHT=137.03598950     TAU= 3.500E+00

 INTERNAL BASIS SET UNCONTRACTION FOR RESC/DK/IOTC/LUT-IOTC:
    44 SHELLS AND     78 AOS ARE CONSTRUCTED, COMPARED TO
    24 SHELLS AND     50 AOS IN THE CONTRACTED BASIS SET


     $CONTRL OPTIONS
     ---------------
 SCFTYP=GVB          RUNTYP=ENERGY       EXETYP=RUN     
 MPLEVL=       0     CITYP =NONE         CCTYP =NONE         VBTYP =NONE    
 DFTTYP=NONE         TDDFT =NONE    
 MULT  =       1     ICHARG=       0     NZVAR =       0     COORD =UNIQUE  
 PP    =NONE         RELWFN=DK           LOCAL =NONE         NUMGRD=       F
 ISPHER=       1     NOSYM =       1     MAXIT =     500     UNITS =ANGS    
 PLTORB=       F     MOLPLT=       F     AIMPAC=       F     FRIEND=        
 NPRINT=       7     IREST =       0     GEOM  =INPUT   
 NORMF =       0     NORMP =       0     ITOL  =      20     ICUT  =      10
 INTTYP=BEST         GRDTYP=BEST         QMTTOL= 1.0E-06

     $SYSTEM OPTIONS
     ---------------
  REPLICATED MEMORY=   125000000 WORDS (ON EVERY NODE).
 DISTRIBUTED MEMDDI=           0 MILLION WORDS IN AGGREGATE,
 MEMDDI DISTRIBUTED OVER   1 PROCESSORS IS           0 WORDS/PROCESSOR.
 TOTAL MEMORY REQUESTED ON EACH PROCESSOR=   125000000 WORDS.
 TIMLIM=      525600.00 MINUTES, OR     365.0 DAYS.
 PARALL= F  BALTYP=  DLB     KDIAG=    0  COREFL= F
 MXSEQ2=     300 MXSEQ3=     150  mem10=         0  mem22=         0

 Using Gaussians for PCM contribution to MEP, params: 0.10000E+01 0.10000E-07

          ----------------
          PROPERTIES INPUT
          ----------------

     MOMENTS            FIELD           POTENTIAL          DENSITY
 IEMOM =       1   IEFLD =       0   IEPOT =       0   IEDEN =       0
 WHERE =COMASS     WHERE =NUCLEI     WHERE =NUCLEI     WHERE =NUCLEI  
 OUTPUT=BOTH       OUTPUT=BOTH       OUTPUT=BOTH       OUTPUT=BOTH    
 IEMINT=       0   IEFINT=       0                     IEDINT=       0
                                                       MORB  =       0
          EXTRAPOLATION IN EFFECT
          SOSCF IN EFFECT
 ORBITAL PRINTING OPTION: NPREO=     1    50     2     1

          *************************
          ROHF-GVB INPUT PARAMETERS
          *************************

          NORB   =   14          NCO    =    2
          NPAIR  =    6          NSETO  =    0
          PAIR ORBITALS
          PAIR    1 HAS ORBS    3    4
          PAIR    2 HAS ORBS    5    6
          PAIR    3 HAS ORBS    7    8
          PAIR    4 HAS ORBS    9   10
          PAIR    5 HAS ORBS   11   12
          PAIR    6 HAS ORBS   13   14

          ----------------------------
          ROHF-GVB COUPLING PARAMETERS
          ----------------------------

          F VECTOR (OCCUPANCIES)
    1   1.0000000000
    2   0.9529411765
    3   0.0470588235
    4   0.9529411765
    5   0.0470588235
    6   0.9529411765
    7   0.0470588235
    8   0.9529411765
    9   0.0470588235
   10   0.9529411765
   11   0.0470588235
   12   0.9529411765
   13   0.0470588235
          ALPHA COUPLING COEFFICEINTS

             1           2           3           4           5

    1    2.0000000
    2    1.9058824   0.9529412
    3    0.0941176   0.0000000   0.0470588
    4    1.9058824   1.8161938   0.0896886   0.9529412
    5    0.0941176   0.0896886   0.0044291   0.0000000   0.0470588
    6    1.9058824   1.8161938   0.0896886   1.8161938   0.0896886
    7    0.0941176   0.0896886   0.0044291   0.0896886   0.0044291
    8    1.9058824   1.8161938   0.0896886   1.8161938   0.0896886
    9    0.0941176   0.0896886   0.0044291   0.0896886   0.0044291
   10    1.9058824   1.8161938   0.0896886   1.8161938   0.0896886
   11    0.0941176   0.0896886   0.0044291   0.0896886   0.0044291
   12    1.9058824   1.8161938   0.0896886   1.8161938   0.0896886
   13    0.0941176   0.0896886   0.0044291   0.0896886   0.0044291

             6           7           8           9          10

    6    0.9529412
    7    0.0000000   0.0470588
    8    1.8161938   0.0896886   0.9529412
    9    0.0896886   0.0044291   0.0000000   0.0470588
   10    1.8161938   0.0896886   1.8161938   0.0896886   0.9529412
   11    0.0896886   0.0044291   0.0896886   0.0044291   0.0000000
   12    1.8161938   0.0896886   1.8161938   0.0896886   1.8161938
   13    0.0896886   0.0044291   0.0896886   0.0044291   0.0896886

            11          12          13

   11    0.0470588
   12    0.0896886   0.9529412
   13    0.0044291   0.0000000   0.0470588
           BETA COUPLING COEFFICIENTS

             1           2           3           4           5

    1   -1.0000000
    2   -0.9529412   0.0000000
    3   -0.0470588  -0.2117647   0.0000000
    4   -0.9529412  -0.9080969  -0.0448443   0.0000000
    5   -0.0470588  -0.0448443  -0.0022145  -0.2117647   0.0000000
    6   -0.9529412  -0.9080969  -0.0448443  -0.9080969  -0.0448443
    7   -0.0470588  -0.0448443  -0.0022145  -0.0448443  -0.0022145
    8   -0.9529412  -0.9080969  -0.0448443  -0.9080969  -0.0448443
    9   -0.0470588  -0.0448443  -0.0022145  -0.0448443  -0.0022145
   10   -0.9529412  -0.9080969  -0.0448443  -0.9080969  -0.0448443
   11   -0.0470588  -0.0448443  -0.0022145  -0.0448443  -0.0022145
   12   -0.9529412  -0.9080969  -0.0448443  -0.9080969  -0.0448443
   13   -0.0470588  -0.0448443  -0.0022145  -0.0448443  -0.0022145

             6           7           8           9          10

    6    0.0000000
    7   -0.2117647   0.0000000
    8   -0.9080969  -0.0448443   0.0000000
    9   -0.0448443  -0.0022145  -0.2117647   0.0000000
   10   -0.9080969  -0.0448443  -0.9080969  -0.0448443   0.0000000
   11   -0.0448443  -0.0022145  -0.0448443  -0.0022145  -0.2117647
   12   -0.9080969  -0.0448443  -0.9080969  -0.0448443  -0.9080969
   13   -0.0448443  -0.0022145  -0.0448443  -0.0022145  -0.0448443

            11          12          13

   11    0.0000000
   12   -0.0448443   0.0000000
   13   -0.0022145  -0.2117647   0.0000000

          NATURAL ORBITAL COEFFICIENTS
          N.O.                PAIR CICOEF-S
            1        0.9761870602    -0.2169304578
            2        0.9761870602    -0.2169304578
            3        0.9761870602    -0.2169304578
            4        0.9761870602    -0.2169304578
            5        0.9761870602    -0.2169304578
            6        0.9761870602    -0.2169304578

     -------------------------------
     INTEGRAL TRANSFORMATION OPTIONS
     -------------------------------
     NWORD  =            0
     CUTOFF = 1.0E-09     MPTRAN =       0
     DIRTRF =       T     AOINTS =DUP     

          ----------------------
          INTEGRAL INPUT OPTIONS
          ----------------------
 NOPK  =       1 NORDER=       0 SCHWRZ=       T

     ------------------------------------------
     THE POINT GROUP IS C1 , NAXIS= 1, ORDER= 1
     ------------------------------------------

 -- VARIATIONAL SPACE WILL BE RESTRICTED TO PURE SPHERICAL HARMONICS ONLY --
 AFTER EXCLUDING CONTAMINANT COMBINATIONS FROM THE CARTESIAN GAUSSIAN BASIS
 SET, THE NUMBER OF SPHERICAL HARMONICS KEPT IN THE VARIATION SPACE IS   48

     DIMENSIONS OF THE SYMMETRY SUBSPACES ARE
 A   =   48

 ..... DONE SETTING UP THE RUN .....
 STEP CPU TIME =     0.02 TOTAL CPU TIME =          0.0 (      0.0 MIN)
 TOTAL WALL CLOCK TIME=          0.1 SECONDS, CPU UTILIZATION IS    13.33%

          ********************
          1 ELECTRON INTEGRALS
          ********************
 TIME TO DO ORDINARY INTEGRALS=      0.01

 ------------------------------------------------------------------
 1ST ORDER DOUGLAS-KROLL-HESS RELATIVISTIC ONE-ELECTRON HAMILTONIAN
            CODED BY TAKAHITO NAKAJIMA AND DMITRI FEDOROV
 ------------------------------------------------------------------

 TIME TO DO  RELATIVISTIC INTS=      0.00
 ...... END OF ONE-ELECTRON INTEGRALS ......
 STEP CPU TIME =     0.01 TOTAL CPU TIME =          0.0 (      0.0 MIN)
 TOTAL WALL CLOCK TIME=          0.2 SECONDS, CPU UTILIZATION IS    12.50%

          -------------
          GUESS OPTIONS
          -------------
          GUESS =MOREAD            NORB  =      48          NORDER=       0
          MIX   =       F          PRTMO =       F          PUNMO =       F
          TOLZ  = 1.0E-08          TOLE  = 1.0E-05
          SYMDEN=       F          PURIFY=       F

 INITIAL GUESS ORBITALS GENERATED BY MOREAD   ROUTINE.

 STATISTICS FOR GENERATION OF SYMMETRY ORBITAL -Q- MATRIX
 NUMBER OF CARTESIAN ATOMIC ORBITALS=         50
 NUMBER OF SPHERICAL CONTAMINANTS DROPPED=     2
 NUMBER OF LINEARLY DEPENDENT MOS DROPPED=     0
 TOTAL NUMBER OF MOS IN VARIATION SPACE=      48

 SYMMETRIES FOR INITIAL GUESS ORBITALS FOLLOW.   BOTH SET(S).
    14 ORBITALS ARE OCCUPIED (    2 CORE ORBITALS).
     1=A        2=A        3=A        4=A        5=A        6=A        7=A   
     8=A        9=A       10=A       11=A       12=A       13=A       14=A   
    15=A       16=A       17=A       18=A       19=A       20=A       21=A   
    22=A       23=A       24=A       25=A       26=A       27=A       28=A   
    29=A       30=A       31=A       32=A       33=A       34=A       35=A   
    36=A       37=A       38=A       39=A       40=A       41=A       42=A   
    43=A       44=A       45=A       46=A       47=A       48=A   
 ...... END OF INITIAL ORBITAL SELECTION ......
 STEP CPU TIME =     0.00 TOTAL CPU TIME =          0.0 (      0.0 MIN)
 TOTAL WALL CLOCK TIME=          0.3 SECONDS, CPU UTILIZATION IS    11.54%

                    ----------------------
                    AO INTEGRAL TECHNOLOGY
                    ----------------------
     S,P,L SHELL ROTATED AXIS INTEGRALS, REPROGRAMMED BY
        KAZUYA ISHIMURA (IMS) AND JOSE SIERRA (SYNSTAR).
     S,P,D,L SHELL ROTATED AXIS INTEGRALS PROGRAMMED BY
        KAZUYA ISHIMURA (INSTITUTE FOR MOLECULAR SCIENCE).
     S,P,D,F,G SHELL TO TOTAL QUARTET ANGULAR MOMENTUM SUM 5,
        ERIC PROGRAM BY GRAHAM FLETCHER (ELORET AND NASA ADVANCED
        SUPERCOMPUTING DIVISION, AMES RESEARCH CENTER).
     S,P,D,F,G,L SHELL GENERAL RYS QUADRATURE PROGRAMMED BY
        MICHEL DUPUIS (PACIFIC NORTHWEST NATIONAL LABORATORY).

          --------------------
          2 ELECTRON INTEGRALS
          --------------------

 DIRECT SCF METHOD SKIPS INTEGRAL STORAGE ON DISK.
 DIRECT TRANSFORMATION SKIPS AO INTEGRAL STORAGE ON DISK.
  ...... END OF TWO-ELECTRON INTEGRALS .....
 STEP CPU TIME =     0.01 TOTAL CPU TIME =          0.0 (      0.0 MIN)
 TOTAL WALL CLOCK TIME=          0.3 SECONDS, CPU UTILIZATION IS    14.81%

          ------------------------
          ROHF-GVB SCF CALCULATION
          ------------------------
 GVB STEP WILL USE    116070 WORDS OF MEMORY.

     MAXIT= 500   NPUNCH= 2   SQCDF TOL=2.0000E-05
     NUCLEAR ENERGY=       31.6688684670
     EXTRAP=T   DAMP=F   SHIFT=F   RSTRCT=F   DIIS=F  SOSCF=T
     SOSCF WILL OPTIMIZE     566 ORBITAL ROTATIONS, SOGTOL=   0.250

 DIRECT SCF CALCULATION, SCHWRZ=T   FDIFF=F
 SCHWARZ INEQUALITY OVERHEAD:      1275 INTEGRALS, T=        0.00

                                                                           NONZERO    BLOCKS
 ITER EX     TOTAL ENERGY       E CHANGE        SQCDF       ORB. GRAD       INTEGRALS   SKIPPED
   0  0      -77.893931794   -77.893931794   1.085111684   0.000000000         750460        55
          ---------------START SECOND ORDER SCF---------------
   1  1      -77.938753639    -0.044821845   0.058415557   0.006487018         750388        63
   2  2      -77.946222473    -0.007468835   0.924667007   0.007109738         750388        63
   3  3      -77.950798733    -0.004576259   0.482463396   0.023156247         750460        55
   4  4      -77.954288115    -0.003489382   0.121942544   0.013089993         750442        57
 SOSCF IS SCALING ROTATION ANGLE MATRIX, SQCDF=  307.025230
   5  5      -77.975834989    -0.021546874  12.684253380   0.015671054         750424        59
 *** RESETTTING SOSCF, UPON ENCOUNTERING A HUGE TOTAL ROTATION=   2.508E+02
   6  6      -62.689988991    15.285845997  12.566911955   1.671892234         750894        27
   7  7      -47.775536187    14.914452804   3.593362719   2.286728946         750889        24
   8  8      -65.097184424   -17.321648236   3.290826517   2.092924835         750898        22
   9  9      -69.904067840    -4.806883416   3.062264869   1.244591263         750895        26
  10 10      -72.881910125    -2.977842285   2.956358187   0.847884143         750894        27
  11 11      -74.103181635    -1.221271510   2.705443669   0.800469069         750893        28
  12 12      -74.898814906    -0.795633272   2.624356622   0.944808661         750894        27
  13 13      -75.289407223    -0.390592316   2.397445466   1.204196617         750895        26
  14 14      -75.564191806    -0.274784584   2.310119948   0.968379700         750894        27
  15 15      -76.154126296    -0.589934490   2.116922900   0.624673646         750815        34
  16 16      -76.301329772    -0.147203476   1.799558615   1.154284918         750813        36
  17 17      -76.605933665    -0.304603893   1.649353621   0.659028898         750754        36
  18 18      -76.521917877     0.084015788   1.566865773   0.832927010         750724        42
  19 19      -75.946612202     0.575305675   1.509682100   1.389025735         750767        35
  20 20      -75.716691097     0.229921106   1.640821842   1.399607080         750711        43
  21 21      -75.634860416     0.081830681   1.753412970   1.869186279         750602        44
  22 22      -76.404279060    -0.769418644   1.544597785   1.030163852         750460        55
  23 23      -76.727814072    -0.323535012   1.182046361   0.691367022         750772        41
  24 24      -77.318016054    -0.590201982   0.791423589   0.563787644         750460        55
  25 25      -77.480872168    -0.162856114   0.715418659   0.489496288         750460        55
  26 26      -77.690623093    -0.209750925   0.549291229   0.284515824         750460        55
  27 27      -77.744750630    -0.054127537   0.541337326   0.316475397         750460        55
          ---------------START SECOND ORDER SCF---------------
  28 28      -77.805805883    -0.061055254   0.212012059   0.191679502         750460        55
  29 29      -77.896876108    -0.091070225   0.057661960   0.074581616         750397        62
  30 30      -77.933190879    -0.036314771   0.077457987   0.036875885         750397        62
  31 31      -77.949997897    -0.016807018   0.992438740   0.011352002         750424        59
  32 32      -77.931751204     0.018246694   0.610919056   0.073165159         750623        47
  33 33      -77.959620723    -0.027869519   0.112186455   0.029492866         750460        55
  34 34      -77.981487216    -0.021866494   0.640992902   0.017316743         750424        59
  35 35      -77.919046438     0.062440778   0.453760589   0.032908529         750583        51
  36 36      -77.987594386    -0.068547948   0.078750654   0.016070918         750460        55
  37 37      -77.989385879    -0.001791493   0.024894600   0.014676361         750388        63
  38 38      -77.997063850    -0.007677971   0.010570229   0.006681426         750388        63
  39 39      -77.998229154    -0.001165304   0.004701939   0.004564163         750442        57
  40 40      -77.998948771    -0.000719617   0.001155762   0.002314636         750397        62
  41 41      -77.999127420    -0.000178649   0.000100909   0.001430531         750388        63
  42 42      -77.999171890    -0.000044471   0.000084888   0.001171214         750388        63
  43 43      -77.999193243    -0.000021352   0.000083583   0.000951679         750388        63
  44 44      -77.999205058    -0.000011815   0.000057478   0.000577729         750388        63
  45 45      -77.999211840    -0.000006782   0.000031064   0.000251206         750388        63
  46 46      -77.999215093    -0.000003253   0.000005168   0.000111966         750388        63
  47 47      -77.999215777    -0.000000684   0.000000867   0.000050314         750388        63
  48 48      -77.999215873    -0.000000096   0.000000154   0.000016160         750388        63
  49 49      -77.999215889    -0.000000016   0.000000065   0.000005807         750388        63
  50 50      -77.999215892    -0.000000004   0.000000022   0.000002497         750388        63
  51 51      -77.999215893    -0.000000001   0.000000007   0.000001496         750388        63

          -----------------
          DENSITY CONVERGED
          -----------------

 FINAL GVB ENERGY IS      -77.9992158929 AFTER  51 ITERATIONS

          ----------------------------
          SCF STATISTICS PER ITERATION
          ----------------------------
          NUMBER OF INTEGRAL PASSES       1
          FOCK FORMATION TIME         0.128
          GEMINAL OPT TIME            0.000
          MIXORB  OPT TIME            0.002
          OCBSE   OPT TIME            0.002

          ----------------
          PAIR INFORMATION
          ----------------
      ORBITAL   CI COEFFICIENTS    OCCUPATION NUMBERS   GVB      ENERGY
 PAIR 1   2     ORB 1     ORB 2      ORB 1     ORB 2    OVERLAP  LOWERING
  1   3   4    0.996348 -0.085382   1.98542   0.01458   0.84214  -0.01424
  2   5   6    0.996348 -0.085383   1.98542   0.01458   0.84214  -0.01424
  3   7   8    0.707125 -0.707088   1.00005   0.99995   0.00003  -0.11916
  4   9  10    0.997325 -0.073100   1.98931   0.01069   0.86342  -0.01242
  5  11  12    0.996348 -0.085383   1.98542   0.01458   0.84214  -0.01424
  6  13  14    0.996348 -0.085383   1.98542   0.01458   0.84214  -0.01424

     THE MAXIMUM LAGRANGIAN ASYMMETRY IS 2.5401821E-06

          ------------
          EIGENVECTORS
          ------------

                      1          2          3          4          5
                  -11.2253   -11.2252    -0.6821    -0.0001    -0.6821
                     A          A          A          A          A   
    1  C  1  S    0.704993  -0.700078   0.014022  -0.004917  -0.136388
    2  C  1  S    0.018907  -0.022729  -0.036667   0.002414   0.279035
    3  C  1  S   -0.014101   0.020456  -0.028067  -0.058613   0.112664
    4  C  1  X    0.001932   0.001872   0.034526  -0.020608   0.180788
    5  C  1  Y    0.000000  -0.000000   0.027830   0.037446  -0.217712
    6  C  1  Z   -0.000000  -0.000000   0.027829   0.037442   0.217725
    7  C  1  X   -0.000749  -0.001617   0.025387   0.020277   0.064380
    8  C  1  Y   -0.000000   0.000000   0.001923   0.053135  -0.089885
    9  C  1  Z   -0.000001   0.000000   0.001923   0.053132   0.089890
   10  C  1 XX   -0.000347  -0.001120   0.001914   0.004379  -0.004635
   11  C  1 YY    0.000174   0.000560  -0.000957  -0.002189   0.002316
   12  C  1 ZZ    0.000174   0.000560  -0.000957  -0.002190   0.002320
   13  C  1 XY    0.000000  -0.000000  -0.006156   0.003581  -0.010122
   14  C  1 XZ    0.000000  -0.000000  -0.006155   0.003581   0.010123
   15  C  1 YZ   -0.000717   0.000863   0.000514   0.001373  -0.009428
   16  H  2  S   -0.002879   0.003207   0.003692  -0.004314   0.425922
   17  H  2  S    0.003802  -0.002877   0.006028  -0.007487   0.128387
   18  H  2  X    0.000234  -0.000639   0.000734   0.001339  -0.013860
   19  H  2  Y   -0.000531   0.000705   0.000517   0.000776   0.015619
   20  H  2  Z    0.000531  -0.000705   0.000104   0.001874  -0.015621
   21  H  3  S   -0.002879   0.003207   0.003690  -0.004306  -0.032251
   22  H  3  S    0.003802  -0.002877   0.006022  -0.007483  -0.034576
   23  H  3  X    0.000234  -0.000639   0.000733   0.001339   0.002605
   24  H  3  Y    0.000531  -0.000705   0.000104   0.001873  -0.006036
   25  H  3  Z   -0.000531   0.000705   0.000517   0.000777   0.006036
   26  C  4  S    0.701178   0.703899  -0.136391   0.184540   0.014024
   27  C  4  S    0.018784   0.022831   0.279041  -0.785524  -0.036672
   28  C  4  S   -0.013989  -0.020533   0.112667   0.169785  -0.028067
   29  C  4  X   -0.001942   0.001861  -0.180762   0.423074  -0.034541
   30  C  4  Y    0.000000  -0.000000   0.217742  -0.507290  -0.027828
   31  C  4  Z   -0.000000  -0.000000   0.217711  -0.507203   0.027830
   32  C  4  X    0.000758  -0.001613  -0.064363   0.060770  -0.025391
   33  C  4  Y   -0.000000  -0.000000   0.089898  -0.059393  -0.001923
   34  C  4  Z    0.000000   0.000000   0.089883  -0.059371   0.001921
   35  C  4 XX   -0.000354   0.001118  -0.004635   0.020076   0.001914
   36  C  4 YY    0.000177  -0.000559   0.002320  -0.010053  -0.000957
   37  C  4 ZZ    0.000177  -0.000559   0.002315  -0.010023  -0.000957
   38  C  4 XY   -0.000000  -0.000000  -0.010124   0.072302  -0.006155
   39  C  4 XZ    0.000000  -0.000000  -0.010121   0.072288   0.006155
   40  C  4 YZ    0.000712   0.000867   0.009428  -0.084458  -0.000514
   41  H  5  S   -0.002861  -0.003223  -0.032267   0.029929   0.003696
   42  H  5  S    0.003786   0.002897  -0.034577  -0.052509   0.006028
   43  H  5  X   -0.000231  -0.000640  -0.002605   0.007974  -0.000734
   44  H  5  Y   -0.000527  -0.000708   0.006036   0.006898  -0.000104
   45  H  5  Z   -0.000527  -0.000708   0.006035   0.006901   0.000517
   46  H  6  S   -0.002861  -0.003223   0.425921   0.980534   0.003695
   47  H  6  S    0.003786   0.002898   0.128387   0.131460   0.006031
   48  H  6  X   -0.000231  -0.000640   0.013861  -0.011150  -0.000734
   49  H  6  Y    0.000527   0.000708  -0.015621   0.016518  -0.000517
   50  H  6  Z    0.000527   0.000708  -0.015619   0.016528   0.000104

                      6          7          8          9         10
                   -0.0001    -0.0906    -0.0906    -0.7672    -0.0001
                     A          A          A          A          A   
    1  C  1  S    0.184537  -0.000016  -0.000010  -0.141423   0.194575
    2  C  1  S   -0.785512   0.000038   0.000024   0.306256  -0.775513
    3  C  1  S    0.169805   0.000050   0.000032   0.078896   0.224041
    4  C  1  X   -0.423127   0.000004   0.000023  -0.354731   0.843149
    5  C  1  Y    0.507197  -0.322765  -0.287295  -0.000010   0.000025
    6  C  1  Z   -0.507255  -0.287303  -0.322778   0.000013  -0.000045
    7  C  1  X   -0.060800   0.000009   0.000023  -0.128199   0.057240
    8  C  1  Y    0.059364  -0.306679  -0.238369  -0.000003   0.000006
    9  C  1  Z   -0.059395  -0.238381  -0.306695   0.000003  -0.000006
   10  C  1 XX    0.020040  -0.000000  -0.000001   0.024713  -0.144210
   11  C  1 YY   -0.010015  -0.000001   0.000001  -0.012357   0.072106
   12  C  1 ZZ   -0.010025   0.000001   0.000000  -0.012357   0.072105
   13  C  1 XY    0.072307  -0.014228   0.015160  -0.000001  -0.000019
   14  C  1 XZ   -0.072315   0.015161  -0.014228  -0.000002   0.000029
   15  C  1 YZ    0.084435   0.000002   0.000003   0.002798  -0.003579
   16  H  2  S    0.980535   0.069244  -0.069266  -0.036379   0.015543
   17  H  2  S    0.131458   0.055386  -0.055404  -0.029252  -0.048066
   18  H  2  X    0.011127  -0.000123   0.000125  -0.005844  -0.001554
   19  H  2  Y   -0.016526  -0.008399  -0.009518   0.000311   0.010165
   20  H  2  Z    0.016536  -0.009518  -0.008398  -0.000311  -0.010168
   21  H  3  S    0.029950  -0.069268   0.069246  -0.036396   0.015518
   22  H  3  S   -0.052471  -0.055407   0.055376  -0.029252  -0.048117
   23  H  3  X   -0.007969   0.000125  -0.000123  -0.005842  -0.001555
   24  H  3  Y   -0.006898  -0.008398  -0.009518  -0.000310  -0.010175
   25  H  3  Z    0.006896  -0.009518  -0.008399   0.000310   0.010174
   26  C  4  S   -0.004911   0.000016  -0.000011  -0.141420  -0.194571
   27  C  4  S    0.002399  -0.000034   0.000022   0.306252   0.775494
   28  C  4  S   -0.058655  -0.000059   0.000054   0.078892  -0.224040
   29  C  4  X    0.020582   0.000007  -0.000021   0.354733   0.843155
   30  C  4  Y   -0.037436   0.322766  -0.287297  -0.000005  -0.000011
   31  C  4  Z    0.037440  -0.287304   0.322781  -0.000006  -0.000027
   32  C  4  X   -0.020313   0.000009  -0.000005   0.128199   0.057245
   33  C  4  Y   -0.053126   0.306678  -0.238365  -0.000001  -0.000004
   34  C  4  Z    0.053133  -0.238379   0.306694   0.000000  -0.000005
   35  C  4 XX    0.004383   0.000001  -0.000004   0.024714   0.144214
   36  C  4 YY   -0.002192  -0.000000   0.000002  -0.012357  -0.072107
   37  C  4 ZZ   -0.002191  -0.000001   0.000001  -0.012357  -0.072106
   38  C  4 XY    0.003580  -0.014228  -0.015160   0.000001  -0.000000
   39  C  4 XZ   -0.003582  -0.015162  -0.014227  -0.000002  -0.000023
   40  C  4 YZ   -0.001375   0.000003  -0.000002  -0.002798  -0.003579
   41  H  5  S   -0.004311  -0.069247  -0.069266  -0.036385  -0.015532
   42  H  5  S   -0.007485  -0.055384  -0.055393  -0.029252   0.048084
   43  H  5  X   -0.001340  -0.000123  -0.000125   0.005843  -0.001554
   44  H  5  Y   -0.001873   0.008399  -0.009519   0.000310  -0.010167
   45  H  5  Z    0.000776  -0.009518   0.008398   0.000310  -0.010170
   46  H  6  S   -0.004307   0.069269   0.069246  -0.036393  -0.015518
   47  H  6  S   -0.007482   0.055410   0.055377  -0.029253   0.048110
   48  H  6  X   -0.001339   0.000124   0.000124   0.005842  -0.001555
   49  H  6  Y   -0.000776   0.008398  -0.009518  -0.000310   0.010174
   50  H  6  Z    0.001872  -0.009519   0.008400  -0.000310   0.010172

                     11         12         13         14         15
                   -0.6821    -0.0001    -0.6821    -0.0001     0.4396
                     A          A          A          A          A   
    1  C  1  S    0.014022  -0.004915  -0.136389   0.184538   0.002450
    2  C  1  S   -0.036668   0.002408   0.279037  -0.785511  -0.168181
    3  C  1  S   -0.028071  -0.058637   0.112666   0.169772   1.773641
    4  C  1  X    0.034535  -0.020595   0.180757  -0.423056  -0.112761
    5  C  1  Y   -0.027828  -0.037437   0.217745  -0.507299   0.000011
    6  C  1  Z   -0.027832  -0.037441  -0.217714   0.507207  -0.000021
    7  C  1  X    0.025392   0.020298   0.064357  -0.060776   0.582467
    8  C  1  Y   -0.001922  -0.053129   0.089900  -0.059394  -0.000033
    9  C  1  Z   -0.001923  -0.053130  -0.089885   0.059373   0.000144
   10  C  1 XX    0.001915   0.004381  -0.004635   0.020081   0.006644
   11  C  1 YY   -0.000957  -0.002191   0.002320  -0.010054  -0.003322
   12  C  1 ZZ   -0.000958  -0.002190   0.002315  -0.010027  -0.003322
   13  C  1 XY    0.006155  -0.003580   0.010123  -0.072300   0.000005
   14  C  1 XZ    0.006156  -0.003581  -0.010121   0.072286   0.000002
   15  C  1 YZ    0.000514   0.001375  -0.009428   0.084460  -0.051965
   16  H  2  S    0.003694  -0.004310  -0.032272   0.029922   0.371311
   17  H  2  S    0.006029  -0.007488  -0.034575  -0.052515  -1.334437
   18  H  2  X    0.000734   0.001340   0.002605  -0.007975   0.020370
   19  H  2  Y   -0.000104  -0.001873   0.006035   0.006898  -0.045567
   20  H  2  Z   -0.000517  -0.000776  -0.006035  -0.006901   0.045570
   21  H  3  S    0.003694  -0.004308   0.425922   0.980527   0.371308
   22  H  3  S    0.006028  -0.007483   0.128390   0.131474  -1.334248
   23  H  3  X    0.000734   0.001339  -0.013861   0.011152   0.020368
   24  H  3  Y   -0.000517  -0.000776  -0.015621   0.016524   0.045563
   25  H  3  Z   -0.000104  -0.001873   0.015618  -0.016526  -0.045561
   26  C  4  S   -0.136389   0.184539   0.014021  -0.004918   0.002449
   27  C  4  S    0.279036  -0.785517  -0.036667   0.002419  -0.168180
   28  C  4  S    0.112665   0.169783  -0.028070  -0.058615   1.773658
   29  C  4  X   -0.180776   0.423098  -0.034524   0.020612   0.112768
   30  C  4  Y   -0.217712   0.507199   0.027831   0.037447   0.000001
   31  C  4  Z   -0.217731   0.507266  -0.027829  -0.037442   0.000019
   32  C  4  X   -0.064372   0.060799  -0.025389  -0.020275  -0.582478
   33  C  4  Y   -0.089886   0.059374   0.001924   0.053135   0.000019
   34  C  4  Z   -0.089894   0.059388  -0.001923  -0.053131  -0.000115
   35  C  4 XX   -0.004636   0.020056   0.001915   0.004379   0.006645
   36  C  4 YY    0.002316  -0.010018  -0.000957  -0.002189  -0.003322
   37  C  4 ZZ    0.002320  -0.010038  -0.000957  -0.002190  -0.003323
   38  C  4 XY    0.010122  -0.072297   0.006156  -0.003582  -0.000006
   39  C  4 XZ    0.010123  -0.072310  -0.006155   0.003580   0.000000
   40  C  4 YZ    0.009429  -0.084443  -0.000514  -0.001374   0.051965
   41  H  5  S    0.425921   0.980533   0.003692  -0.004310   0.371301
   42  H  5  S    0.128390   0.131467   0.006028  -0.007488  -1.334380
   43  H  5  X    0.013860  -0.011134  -0.000734  -0.001339  -0.020369
   44  H  5  Y    0.015619  -0.016524   0.000517   0.000777  -0.045566
   45  H  5  Z    0.015621  -0.016533  -0.000104  -0.001873  -0.045568
   46  H  6  S   -0.032257   0.029947   0.003690  -0.004306   0.371326
   47  H  6  S   -0.034574  -0.052477   0.006021  -0.007484  -1.334338
   48  H  6  X   -0.002605   0.007971  -0.000733  -0.001339  -0.020371
   49  H  6  Y   -0.006036  -0.006899   0.000104   0.001873   0.045565
   50  H  6  Z   -0.006036  -0.006897  -0.000517  -0.000777   0.045564

                     16         17         18         19         20
                    0.4745     0.4816     0.4816     0.6047     0.6050
                     A          A          A          A          A   
    1  C  1  S    0.043080  -0.000002  -0.000011   0.103160   0.024473
    2  C  1  S   -0.531849   0.000002   0.000127  -0.845917  -0.107523
    3  C  1  S    3.276887  -0.000209  -0.000729   0.860358  -2.492939
    4  C  1  X   -0.050028  -0.000017   0.000051   0.405851  -0.332125
    5  C  1  Y    0.000032  -0.127662   0.139622  -0.000014   0.000055
    6  C  1  Z   -0.000046  -0.085757  -0.168627   0.000024  -0.000017
    7  C  1  X   -0.050746   0.000245  -0.000204  -0.641806   3.111674
    8  C  1  Y   -0.000217   0.463287  -1.179170   0.000049  -0.000094
    9  C  1  Z    0.000395   0.131652   1.260013  -0.000085   0.000045
   10  C  1 XX   -0.001011  -0.000003  -0.000006  -0.069969  -0.023008
   11  C  1 YY    0.000506   0.000004  -0.000008   0.034984   0.011503
   12  C  1 ZZ    0.000505  -0.000001   0.000013   0.034984   0.011505
   13  C  1 XY   -0.000011  -0.008254  -0.058944   0.000001   0.000002
   14  C  1 XZ    0.000017  -0.023690   0.054600   0.000002   0.000004
   15  C  1 YZ   -0.043272   0.000001   0.000009   0.003464   0.022403
   16  H  2  S    0.310301  -0.057055   0.419551   0.080428  -0.184407
   17  H  2  S   -1.447900   0.296796  -2.183045  -0.282266  -0.565394
   18  H  2  X    0.028390  -0.010638   0.078235  -0.037389   0.010241
   19  H  2  Y   -0.042625  -0.005656  -0.054054  -0.004961  -0.015922
   20  H  2  Z    0.042628  -0.019881   0.050597   0.004963   0.015931
   21  H  3  S    0.310063   0.057044  -0.419746   0.080406  -0.184450
   22  H  3  S   -1.446762  -0.296927   2.183955  -0.282351  -0.565151
   23  H  3  X    0.028354   0.010638  -0.078259  -0.037384   0.010230
   24  H  3  Y    0.042599  -0.005646  -0.054096   0.004967   0.015917
   25  H  3  Z   -0.042595  -0.019884   0.050611  -0.004965  -0.015907
   26  C  4  S   -0.043081   0.000001   0.000009   0.103203  -0.024291
   27  C  4  S    0.531852  -0.000016  -0.000135  -0.846102   0.106039
   28  C  4  S   -3.276876   0.000239   0.000662   0.855929   2.494418
   29  C  4  X   -0.050026  -0.000032   0.000021  -0.405266  -0.332843
   30  C  4  Y   -0.000013   0.168628  -0.085756  -0.000014  -0.000048
   31  C  4  Z    0.000027   0.139616   0.127663  -0.000024  -0.000014
   32  C  4  X   -0.050751   0.000262  -0.000151   0.636301   3.112818
   33  C  4  Y   -0.000066  -1.260045   0.131670   0.000052   0.000089
   34  C  4  Z   -0.000062  -1.179132  -0.463266   0.000082   0.000035
   35  C  4 XX    0.001010  -0.000002  -0.000004  -0.070009   0.022884
   36  C  4 YY   -0.000506  -0.000008   0.000002   0.035004  -0.011441
   37  C  4 ZZ   -0.000505   0.000011   0.000003   0.035005  -0.011442
   38  C  4 XY    0.000010   0.054601   0.023691  -0.000001  -0.000002
   39  C  4 XZ   -0.000005   0.058943  -0.008254   0.000002  -0.000004
   40  C  4 YZ   -0.043271   0.000004   0.000009  -0.003503   0.022398
   41  H  5  S   -0.310172   0.419656   0.057116   0.080102   0.184558
   42  H  5  S    1.447227  -2.183456  -0.297242  -0.283263   0.564902
   43  H  5  X    0.028368  -0.078245  -0.010644   0.037370   0.010307
   44  H  5  Y    0.042604  -0.050595  -0.019895  -0.004989   0.015912
   45  H  5  Z    0.042615  -0.054080   0.005639  -0.004991   0.015921
   46  H  6  S   -0.310183  -0.419642  -0.056993   0.080081   0.184602
   47  H  6  S    1.447413   2.183541   0.296456  -0.283344   0.564667
   48  H  6  X    0.028377   0.078250   0.010631   0.037366   0.010297
   49  H  6  Y   -0.042617  -0.050611  -0.019870   0.004995  -0.015908
   50  H  6  Z   -0.042607  -0.054070   0.005663   0.004993  -0.015899

                     21         22         23         24         25
                    0.6957     0.7267     0.7527     0.7527     0.8028
                     A          A          A          A          A   
    1  C  1  S    0.103791  -0.092313  -0.000003  -0.000004   0.000016
    2  C  1  S   -0.887205   1.048522   0.000035  -0.000038  -0.000154
    3  C  1  S    2.754421  -5.637938  -0.000225   0.000492   0.000363
    4  C  1  X   -0.216968  -0.000244  -0.000002   0.000014   0.000013
    5  C  1  Y   -0.000004  -0.000052  -0.229516  -0.660097   0.145888
    6  C  1  Z    0.000042  -0.000043   0.042115  -0.697591  -0.331496
    7  C  1  X    1.085137   1.925859   0.000119  -0.000107   0.000062
    8  C  1  Y    0.000120   0.000072   0.504797   0.627273  -0.822861
    9  C  1  Z   -0.000082  -0.000077  -0.316010   0.740489   1.784693
   10  C  1 XX    0.006535  -0.012528  -0.000001   0.000001  -0.000005
   11  C  1 YY   -0.003269   0.006265  -0.000000  -0.000001   0.000016
   12  C  1 ZZ   -0.003266   0.006263   0.000002   0.000000  -0.000011
   13  C  1 XY   -0.000011  -0.000001   0.035221   0.007724   0.082695
   14  C  1 XZ    0.000000   0.000012  -0.031813   0.016983  -0.007291
   15  C  1 YZ    0.028659  -0.058340  -0.000005   0.000005   0.000009
   16  H  2  S   -0.347061   0.395269   0.042394  -0.005873  -0.482784
   17  H  2  S   -1.013729   0.292021   0.616130  -0.085117  -0.627070
   18  H  2  X    0.025999  -0.024766   0.011775  -0.001616  -0.007023
   19  H  2  Y   -0.011692  -0.008165   0.026203  -0.030417  -0.012775
   20  H  2  Z    0.011700   0.008156  -0.033466  -0.022175   0.048288
   21  H  3  S   -0.347165   0.395168  -0.042358   0.005772   0.482652
   22  H  3  S   -1.013822   0.292075  -0.616152   0.084948   0.626859
   23  H  3  X    0.025999  -0.024766  -0.011772   0.001632   0.007031
   24  H  3  Y    0.011704   0.008153   0.026209  -0.030411  -0.012768
   25  H  3  Z   -0.011697  -0.008163  -0.033467  -0.022185   0.048308
   26  C  4  S    0.103789   0.092316  -0.000006  -0.000003   0.000015
   27  C  4  S   -0.887180  -1.048546  -0.000013   0.000051  -0.000152
   28  C  4  S    2.754339   5.638000   0.000361  -0.000385   0.000330
   29  C  4  X    0.216960  -0.000239  -0.000017   0.000002  -0.000002
   30  C  4  Y   -0.000002   0.000052  -0.697591  -0.042107   0.122079
   31  C  4  Z   -0.000042  -0.000037   0.660092  -0.229522   0.321720
   32  C  4  X   -1.085121   1.925838   0.000063  -0.000148  -0.000106
   33  C  4  Y    0.000118  -0.000070   0.740559   0.315978  -0.694738
   34  C  4  Z    0.000079  -0.000089  -0.627197   0.504800  -1.729354
   35  C  4 XX    0.006536   0.012528   0.000001  -0.000002  -0.000005
   36  C  4 YY   -0.003270  -0.006264  -0.000001   0.000000   0.000016
   37  C  4 ZZ   -0.003267  -0.006263   0.000001   0.000002  -0.000010
   38  C  4 XY    0.000012  -0.000003  -0.016978  -0.031810  -0.082383
   39  C  4 XZ    0.000001  -0.000012   0.007728  -0.035226  -0.001272
   40  C  4 YZ   -0.028658  -0.058340  -0.000004   0.000006  -0.000010
   41  H  5  S   -0.347047  -0.395273   0.005815   0.042399  -0.448822
   42  H  5  S   -1.013716  -0.292046   0.085001   0.616151  -0.582967
   43  H  5  X   -0.025999  -0.024768  -0.001633  -0.011772   0.006527
   44  H  5  Y   -0.011692   0.008164  -0.022180   0.033463  -0.009282
   45  H  5  Z   -0.011700   0.008155   0.030413   0.026208  -0.047485
   46  H  6  S   -0.347160  -0.395174  -0.005906  -0.042338   0.448685
   47  H  6  S   -1.013793  -0.292098  -0.085145  -0.616113   0.582737
   48  H  6  X   -0.025998  -0.024766   0.001617   0.011775  -0.006537
   49  H  6  Y    0.011703  -0.008153  -0.022174   0.033467  -0.009275
   50  H  6  Z    0.011696  -0.008163   0.030422   0.026205  -0.047505

                     26         27         28         29         30
                    0.8028     1.4434     1.4858     1.7667     1.7668
                     A          A          A          A          A   
    1  C  1  S   -0.000007   0.000001   0.000001   0.000018  -0.000002
    2  C  1  S    0.000086  -0.000009  -0.000002  -0.000127   0.000028
    3  C  1  S   -0.000244   0.000019  -0.000012  -0.000109  -0.000085
    4  C  1  X    0.000016  -0.000005  -0.000004  -0.000016   0.000005
    5  C  1  Y    0.321696  -0.000000   0.000007  -0.020534  -0.026371
    6  C  1  Z   -0.122103   0.000006  -0.000010   0.025464   0.021643
    7  C  1  X   -0.000098   0.000009   0.000007   0.000412  -0.000058
    8  C  1  Y   -1.729374   0.000006  -0.000059   0.408556  -0.024317
    9  C  1  Z    0.694763  -0.000011   0.000035   0.041758  -0.407143
   10  C  1 XX   -0.000003  -0.000000  -0.000001  -0.000035  -0.000004
   11  C  1 YY    0.000007   0.278850   0.235037   0.000027   0.000011
   12  C  1 ZZ   -0.000004  -0.278850  -0.235036   0.000008  -0.000008
   13  C  1 XY    0.001272  -0.000068  -0.000008  -0.183009   0.275002
   14  C  1 XZ   -0.082383   0.000005   0.000067  -0.282573   0.171089
   15  C  1 YZ   -0.000011   0.000004   0.000006  -0.000114   0.000029
   16  H  2  S   -0.448704  -0.000002  -0.000013   0.086907   0.090788
   17  H  2  S   -0.582802   0.000010  -0.000028   0.048215   0.050462
   18  H  2  X   -0.006544   0.000050   0.000043   0.260263   0.271759
   19  H  2  Y   -0.047504  -0.254364  -0.277979  -0.108307   0.323436
   20  H  2  Z    0.009279  -0.254376  -0.277981  -0.327736   0.094327
   21  H  3  S    0.448822   0.000006   0.000020  -0.087028  -0.090741
   22  H  3  S    0.582942  -0.000028   0.000027  -0.048424  -0.050394
   23  H  3  X    0.006523  -0.000051  -0.000043  -0.260459  -0.271715
   24  H  3  Y   -0.047490   0.254328   0.278012  -0.108221   0.323419
   25  H  3  Z    0.009273   0.254323   0.278018  -0.327783   0.094376
   26  C  4  S    0.000005   0.000002  -0.000001  -0.000018  -0.000004
   27  C  4  S   -0.000071  -0.000009   0.000001   0.000124   0.000040
   28  C  4  S    0.000210   0.000021   0.000012   0.000121  -0.000118
   29  C  4  X    0.000016   0.000004  -0.000004  -0.000011   0.000004
   30  C  4  Y   -0.331477  -0.000005  -0.000003   0.021642  -0.025470
   31  C  4  Z   -0.145918  -0.000009  -0.000006   0.026367  -0.020536
   32  C  4  X   -0.000090  -0.000007   0.000007   0.000402   0.000046
   33  C  4  Y    1.784721   0.000012   0.000054  -0.407142  -0.041757
   34  C  4  Z    0.822891   0.000015   0.000028   0.024319   0.408557
   35  C  4 XX    0.000003  -0.000000   0.000001   0.000032  -0.000006
   36  C  4 YY   -0.000008   0.278852  -0.235035  -0.000026   0.000012
   37  C  4 ZZ    0.000005  -0.278851   0.235034  -0.000006  -0.000005
   38  C  4 XY    0.007290   0.000067  -0.000006  -0.171089  -0.282573
   39  C  4 XZ    0.082695   0.000003  -0.000067   0.275003   0.183009
   40  C  4 YZ   -0.000010  -0.000004   0.000006  -0.000110  -0.000032
   41  H  5  S    0.482670  -0.000004   0.000014  -0.090706   0.086989
   42  H  5  S    0.626936   0.000015   0.000023  -0.050331   0.048354
   43  H  5  X   -0.007037  -0.000050   0.000042   0.271643  -0.260383
   44  H  5  Y    0.048309  -0.254365   0.277977   0.094383   0.327769
   45  H  5  Z    0.012772   0.254378  -0.277979  -0.323406  -0.108250
   46  H  6  S   -0.482777   0.000008  -0.000021   0.090826  -0.086938
   47  H  6  S   -0.627054  -0.000033  -0.000022   0.050530  -0.048281
   48  H  6  X    0.007017   0.000050  -0.000042  -0.271832   0.260339
   49  H  6  Y    0.048294   0.254330  -0.278011   0.094301   0.327750
   50  H  6  Z    0.012767  -0.254324   0.278017  -0.323449  -0.108297

                     31         32         33         34         35
                    1.8612     1.8733     1.8733     1.9023     2.1338
                     A          A          A          A          A   
    1  C  1  S    0.025189   0.000020  -0.000032   0.057950   0.000072
    2  C  1  S   -0.195192  -0.000131   0.000232  -0.374378  -0.000520
    3  C  1  S   -0.384704  -0.000331   0.000299  -0.934822   0.000286
    4  C  1  X    0.270557   0.000000  -0.000017  -0.049998  -0.000066
    5  C  1  Y    0.000022  -0.216793  -0.017201   0.000065  -0.055525
    6  C  1  Z   -0.000011   0.216950  -0.015000  -0.000110   0.260928
    7  C  1  X   -0.750626   0.000541  -0.000530   1.422308   0.000118
    8  C  1  Y   -0.000059   1.095451  -0.371702  -0.000588   0.052604
    9  C  1  Z   -0.000074  -1.091600  -0.382773   0.000206   0.141683
   10  C  1 XX    0.297311  -0.000044   0.000015  -0.109987   0.000242
   11  C  1 YY   -0.148658   0.000049  -0.000001   0.054999  -0.000091
   12  C  1 ZZ   -0.148653  -0.000005  -0.000014   0.054989  -0.000151
   13  C  1 XY   -0.000005   0.208477   0.201983   0.000039   0.412498
   14  C  1 XZ    0.000035  -0.210512   0.199860   0.000275  -0.063657
   15  C  1 YZ    0.098668  -0.000136   0.000169  -0.383067  -0.000373
   16  H  2  S    0.149708   0.413553   0.002168  -0.167272   0.038605
   17  H  2  S    0.161885   0.496165   0.002618  -0.199946  -0.127549
   18  H  2  X    0.386644  -0.350943  -0.001663  -0.395437   0.256210
   19  H  2  Y    0.103750  -0.314115   0.181365  -0.037710  -0.297250
   20  H  2  Z   -0.103747   0.312112   0.184511   0.037974  -0.182952
   21  H  3  S    0.149690  -0.413706  -0.002015  -0.166911  -0.038728
   22  H  3  S    0.161896  -0.496323  -0.002405  -0.199643   0.127476
   23  H  3  X    0.386651   0.350646   0.001893  -0.395504  -0.256478
   24  H  3  Y   -0.103761  -0.313956   0.181339   0.038180  -0.297532
   25  H  3  Z    0.103765   0.312233   0.184552  -0.037912  -0.182765
   26  C  4  S    0.025185  -0.000029   0.000018  -0.057953   0.000070
   27  C  4  S   -0.195166   0.000205  -0.000115   0.374396  -0.000530
   28  C  4  S   -0.384642   0.000207  -0.000318   0.934839   0.000501
   29  C  4  X   -0.270564  -0.000006   0.000047  -0.049985   0.000095
   30  C  4  Y   -0.000009  -0.015002  -0.216953  -0.000078   0.023715
   31  C  4  Z   -0.000023   0.017197  -0.216784  -0.000120  -0.256054
   32  C  4  X    0.750721   0.000573  -0.000687   1.422263   0.000006
   33  C  4  Y    0.000120  -0.382786   1.091623   0.000610   0.099611
   34  C  4  Z    0.000166   0.371687   1.095423   0.000259  -0.165014
   35  C  4 XX    0.297320   0.000046  -0.000085   0.109969   0.000256
   36  C  4 YY   -0.148662  -0.000017   0.000069  -0.054989  -0.000100
   37  C  4 ZZ   -0.148658  -0.000029   0.000016  -0.054980  -0.000156
   38  C  4 XY    0.000001  -0.199860  -0.210510   0.000026  -0.412361
   39  C  4 XZ   -0.000016   0.201983  -0.208475  -0.000274   0.062767
   40  C  4 YZ   -0.098692  -0.000171   0.000159  -0.383062   0.000377
   41  H  5  S    0.149761  -0.002009   0.413536   0.167269   0.028327
   42  H  5  S    0.161965  -0.002398   0.496132   0.199957  -0.093660
   43  H  5  X   -0.386621  -0.001931   0.351001  -0.395400  -0.188100
   44  H  5  Y    0.103697   0.184550  -0.312276   0.037707  -0.368956
   45  H  5  Z    0.103732  -0.181325  -0.313985   0.037950   0.285111
   46  H  6  S    0.149652   0.002193  -0.413726   0.166891  -0.028457
   47  H  6  S    0.161846   0.002639  -0.496351   0.199617   0.093583
   48  H  6  X   -0.386721   0.001626  -0.350587  -0.395490   0.188319
   49  H  6  Y   -0.103822   0.184500  -0.312087  -0.038172  -0.369248
   50  H  6  Z   -0.103787  -0.181386  -0.314072  -0.037923   0.284911

                     36         37         38         39         40
                    2.1338     2.1892     2.3511     2.4142     2.6107
                     A          A          A          A          A   
    1  C  1  S   -0.000015  -0.083611   0.028255  -0.000001   0.000011
    2  C  1  S    0.000174   0.604582  -0.419773   0.000012  -0.000047
    3  C  1  S   -0.000937  -0.364166   3.066071  -0.000006  -0.000451
    4  C  1  X    0.000074   0.082883  -0.298158   0.000003  -0.000021
    5  C  1  Y   -0.256057  -0.000002  -0.000084  -0.000098  -0.346344
    6  C  1  Z   -0.023720   0.000254   0.000012   0.000082   0.332849
    7  C  1  X    0.000474  -0.034397  -1.903869   0.000003   0.000248
    8  C  1  Y   -0.165006   0.000042  -0.000129  -0.000031  -0.315734
    9  C  1  Z   -0.099617   0.000168   0.000106  -0.000012   0.265414
   10  C  1 XX   -0.000114  -0.338975   0.369798  -0.000008  -0.000111
   11  C  1 YY    0.000059   0.169471  -0.184895   0.591676   0.000735
   12  C  1 ZZ    0.000055   0.169504  -0.184903  -0.591668  -0.000624
   13  C  1 XY   -0.062768   0.000318   0.000040  -0.000099   0.457259
   14  C  1 XZ   -0.412365  -0.000021  -0.000067  -0.000099   0.340165
   15  C  1 YZ    0.000094   0.450368  -0.401531   0.000019   0.000070
   16  H  2  S    0.028417   0.039982  -0.023885  -0.000010  -0.203789
   17  H  2  S   -0.093602   0.009817   0.096257  -0.000056  -0.325316
   18  H  2  X    0.188149   0.177102   0.234047   0.000054   0.363994
   19  H  2  Y    0.284957  -0.162689   0.272326   0.292735  -0.330245
   20  H  2  Z    0.369155   0.162180  -0.272257   0.292835   0.064987
   21  H  3  S   -0.028359   0.039958  -0.023825   0.000008   0.203842
   22  H  3  S    0.093630   0.010040   0.096441   0.000051   0.325355
   23  H  3  X   -0.188259   0.176714   0.233877  -0.000040  -0.364093
   24  H  3  Y    0.285054   0.162094  -0.272313  -0.292709  -0.330921
   25  H  3  Z    0.369060  -0.162574   0.272379  -0.292639   0.064363
   26  C  4  S   -0.000007  -0.083612  -0.028254  -0.000002  -0.000009
   27  C  4  S   -0.000013   0.604585   0.419767   0.000011   0.000031
   28  C  4  S    0.000812  -0.364177  -3.066067   0.000003   0.000382
   29  C  4  X    0.000048  -0.082880  -0.298156  -0.000002   0.000029
   30  C  4  Y    0.260933  -0.000000   0.000088  -0.000099   0.375918
   31  C  4  Z    0.055518  -0.000255   0.000018  -0.000083   0.363521
   32  C  4  X    0.000495   0.034389  -1.903869   0.000004   0.000278
   33  C  4  Y    0.141676   0.000043   0.000126  -0.000031   0.339461
   34  C  4  Z   -0.052610  -0.000168   0.000103   0.000012   0.293240
   35  C  4 XX    0.000039  -0.338976  -0.369796  -0.000007   0.000036
   36  C  4 YY   -0.000030   0.169472   0.184897   0.591676  -0.000681
   37  C  4 ZZ   -0.000008   0.169504   0.184899  -0.591669   0.000645
   38  C  4 XY    0.063656  -0.000318   0.000049   0.000098   0.429965
   39  C  4 XZ    0.412502  -0.000020   0.000064  -0.000099  -0.302486
   40  C  4 YZ   -0.000019  -0.450369  -0.401529  -0.000017   0.000085
   41  H  5  S   -0.038675   0.039982   0.023887  -0.000012   0.221885
   42  H  5  S    0.127509   0.009817  -0.096258  -0.000054   0.354209
   43  H  5  X    0.256246  -0.177100   0.234049  -0.000053   0.396255
   44  H  5  Y   -0.182847  -0.162691  -0.272321   0.292738   0.336992
   45  H  5  Z    0.297419  -0.162181  -0.272261  -0.292833   0.093392
   46  H  6  S    0.038655   0.039959   0.023825   0.000007  -0.221904
   47  H  6  S   -0.127509   0.010039  -0.096442   0.000053  -0.354182
   48  H  6  X   -0.256432  -0.176715   0.233874   0.000040  -0.396418
   49  H  6  Y   -0.182859   0.162093   0.272316  -0.292711   0.337571
   50  H  6  Z    0.297378   0.162576   0.272378   0.292637   0.092688

                     41         42         43         44         45
                    2.6107     2.6351     2.6486     2.9723     3.1246
                     A          A          A          A          A   
    1  C  1  S    0.000041  -0.022755   0.000002  -0.132833  -0.054251
    2  C  1  S   -0.000264   0.152420  -0.000013   0.755012   0.408477
    3  C  1  S   -0.000968   0.497832   0.000012   0.876984   0.632197
    4  C  1  X   -0.000736   0.357381  -0.000008   0.195204   0.477459
    5  C  1  Y    0.363524   0.000756   0.000366  -0.000050   0.000035
    6  C  1  Z   -0.375905  -0.000750  -0.000337   0.000127   0.000007
    7  C  1  X   -0.000482   0.251296  -0.000013   0.082191  -0.115307
    8  C  1  Y    0.293242   0.000584   0.000408  -0.000059   0.000008
    9  C  1  Z   -0.339451  -0.000635  -0.000325   0.000082   0.000020
   10  C  1 XX   -0.001082   0.528182   0.000002   0.658219  -0.718499
   11  C  1 YY    0.000704  -0.264089   0.649477  -0.329099   0.359266
   12  C  1 ZZ    0.000378  -0.264094  -0.649479  -0.329120   0.359233
   13  C  1 XY    0.302487   0.000705  -0.000548  -0.000069   0.000104
   14  C  1 XZ    0.429979   0.000888  -0.000322   0.000108   0.000013
   15  C  1 YZ   -0.000239   0.126610   0.000016  -0.380203  -0.615761
   16  H  2  S    0.222156  -0.138998   0.000247  -0.402625  -0.424426
   17  H  2  S    0.354669  -0.241585   0.000376  -0.288253  -0.114947
   18  H  2  X   -0.395836  -0.240501  -0.000435   0.112887   0.534884
   19  H  2  Y    0.093721  -0.326109   0.275041  -0.432385  -0.209319
   20  H  2  Z   -0.337885   0.325601   0.274622   0.432447   0.209360
   21  H  3  S   -0.221631  -0.139936  -0.000237  -0.402444  -0.424460
   22  H  3  S   -0.353711  -0.242955  -0.000368  -0.288131  -0.114936
   23  H  3  X    0.396825  -0.238942   0.000430   0.112759   0.534910
   24  H  3  Y    0.092359   0.326490  -0.274362   0.432296   0.209398
   25  H  3  Z   -0.336670  -0.326991  -0.274765  -0.432221  -0.209357
   26  C  4  S    0.000041  -0.022755  -0.000001   0.132833   0.054255
   27  C  4  S   -0.000266   0.152420   0.000011  -0.755014  -0.408503
   28  C  4  S   -0.001010   0.497834  -0.000013  -0.876982  -0.632230
   29  C  4  X    0.000735  -0.357382  -0.000006   0.195204   0.477475
   30  C  4  Y    0.332850   0.000742  -0.000367   0.000053  -0.000046
   31  C  4  Z    0.346330   0.000737  -0.000338   0.000131  -0.000004
   32  C  4  X    0.000454  -0.251296  -0.000012   0.082190  -0.115296
   33  C  4  Y    0.265417   0.000573  -0.000410   0.000058  -0.000010
   34  C  4  Z    0.315725   0.000627  -0.000327   0.000081   0.000020
   35  C  4 XX   -0.001087   0.528186  -0.000006  -0.658217   0.718503
   36  C  4 YY    0.000761  -0.264094  -0.649475   0.329097  -0.359268
   37  C  4 ZZ    0.000327  -0.264092   0.649480   0.329120  -0.359235
   38  C  4 XY   -0.340167  -0.000724  -0.000547  -0.000065   0.000105
   39  C  4 XZ    0.457272   0.000899   0.000320  -0.000104  -0.000006
   40  C  4 YZ    0.000232  -0.126610   0.000017  -0.380203  -0.615795
   41  H  5  S    0.204077  -0.139008  -0.000246   0.402628   0.424439
   42  H  5  S    0.325809  -0.241594  -0.000378   0.288252   0.114957
   43  H  5  X    0.363550   0.240491  -0.000440   0.112888   0.534900
   44  H  5  Y    0.065382  -0.326120  -0.275040   0.432382   0.209331
   45  H  5  Z    0.331159  -0.325599   0.274621   0.432445   0.209374
   46  H  6  S   -0.203550  -0.139929   0.000239   0.402443   0.424491
   47  H  6  S   -0.324852  -0.242946   0.000371   0.288131   0.114953
   48  H  6  X   -0.364525   0.238955   0.000430   0.112759   0.534934
   49  H  6  Y    0.063971   0.326483   0.274359  -0.432296  -0.209420
   50  H  6  Z    0.330000   0.326991  -0.274768  -0.432223  -0.209378

                     46         47         48
                    3.2018     3.2677     3.2677
                     A          A          A   
    1  C  1  S    0.118231   0.000008  -0.000018
    2  C  1  S   -0.783888  -0.000058   0.000103
    3  C  1  S   -0.643932   0.000003   0.000037
    4  C  1  X   -0.425043  -0.000046   0.000018
    5  C  1  Y   -0.000005  -0.500706   0.390761
    6  C  1  Z   -0.000076   0.436896  -0.460984
    7  C  1  X   -0.254944   0.000015   0.000090
    8  C  1  Y    0.000001  -0.288022   0.359814
    9  C  1  Z   -0.000034   0.385725  -0.252243
   10  C  1 XX    0.174933   0.000024   0.000168
   11  C  1 YY   -0.087477  -0.000078   0.000023
   12  C  1 ZZ   -0.087456   0.000053  -0.000191
   13  C  1 XY    0.000006  -0.529008   0.854693
   14  C  1 XZ   -0.000076   0.901398  -0.444741
   15  C  1 YZ    0.826041   0.000055  -0.000010
   16  H  2  S    0.590947  -0.624949   0.567702
   17  H  2  S    0.276748  -0.252263   0.229116
   18  H  2  X   -0.477475   0.290964  -0.264395
   19  H  2  Y    0.400008  -0.491936   0.360378
   20  H  2  Z   -0.400095   0.405936  -0.455168
   21  H  3  S    0.590869   0.625019  -0.567815
   22  H  3  S    0.276723   0.252262  -0.229207
   23  H  3  X   -0.477460  -0.291084   0.264349
   24  H  3  Y   -0.400041  -0.492037   0.360544
   25  H  3  Z    0.399931   0.405859  -0.455249
   26  C  4  S    0.118228   0.000006   0.000017
   27  C  4  S   -0.783866  -0.000047  -0.000096
   28  C  4  S   -0.643909   0.000011  -0.000032
   29  C  4  X    0.425021   0.000028   0.000012
   30  C  4  Y   -0.000013  -0.460996  -0.436902
   31  C  4  Z    0.000068  -0.390755  -0.500694
   32  C  4  X    0.254946  -0.000022   0.000088
   33  C  4  Y   -0.000004  -0.252250  -0.385734
   34  C  4  Z    0.000032  -0.359805  -0.288014
   35  C  4 XX    0.174902   0.000035  -0.000164
   36  C  4 YY   -0.087462  -0.000072  -0.000031
   37  C  4 ZZ   -0.087440   0.000037   0.000194
   38  C  4 XY    0.000002   0.444753   0.901413
   39  C  4 XZ   -0.000070   0.854680   0.528996
   40  C  4 YZ   -0.826015  -0.000049  -0.000007
   41  H  5  S    0.590918  -0.567734  -0.624934
   42  H  5  S    0.276740  -0.229169  -0.252217
   43  H  5  X    0.477450  -0.264319  -0.291040
   44  H  5  Y    0.399993  -0.455176  -0.405823
   45  H  5  Z    0.400081  -0.360498  -0.491945
   46  H  6  S    0.590856   0.567780   0.625036
   47  H  6  S    0.276721   0.229159   0.252303
   48  H  6  X    0.477442   0.264434   0.290997
   49  H  6  Y   -0.400035  -0.455256  -0.405988
   50  H  6  Z   -0.399927  -0.360407  -0.492008

          -----------
          GI ORBITALS
          -----------


                    PAIR   1

                      1          2

    1  C  1  S    0.012075  -0.014838
    2  C  1  S   -0.034512   0.035869
    3  C  1  S   -0.043404   0.010470
    4  C  1  X    0.027346  -0.038925
    5  C  1  Y    0.037230  -0.016189
    6  C  1  Z    0.037227  -0.016188
    7  C  1  X    0.030061  -0.018667
    8  C  1  Y    0.016774   0.013083
    9  C  1  Z    0.016773   0.013082
   10  C  1 XX    0.003068  -0.000607
   11  C  1 YY   -0.001534   0.000304
   12  C  1 ZZ   -0.001534   0.000303
   13  C  1 XY   -0.004902   0.006914
   14  C  1 XZ   -0.004901   0.006913
   15  C  1 YZ    0.000879  -0.000108
   16  H  2  S    0.002332  -0.004755
   17  H  2  S    0.003682  -0.007889
   18  H  2  X    0.001080  -0.000328
   19  H  2  Y    0.000714  -0.000278
   20  H  2  Z    0.000626   0.000427
   21  H  3  S    0.002331  -0.004751
   22  H  3  S    0.003677  -0.007881
   23  H  3  X    0.001080  -0.000328
   24  H  3  Y    0.000626   0.000427
   25  H  3  Z    0.000714  -0.000278
   26  C  4  S   -0.079051   0.182743
   27  C  4  S    0.047111  -0.488492
   28  C  4  S    0.155829  -0.060428
   29  C  4  X   -0.054621   0.292343
   30  C  4  Y    0.066450  -0.351493
   31  C  4  Z    0.066445  -0.351439
   32  C  4  X   -0.044698   0.078844
   33  C  4  Y    0.069591  -0.102964
   34  C  4  Z    0.069583  -0.102943
   35  C  4 XX    0.001192   0.010088
   36  C  4 YY   -0.000598  -0.005051
   37  C  4 ZZ   -0.000595  -0.005037
   38  C  4 XY    0.010597   0.030029
   39  C  4 XZ    0.010595   0.030023
   40  C  4 YZ   -0.014680  -0.032777
   41  H  5  S   -0.022558   0.039376
   42  H  5  S   -0.047936   0.018432
   43  H  5  X   -0.000260   0.004740
   44  H  5  Y    0.007730  -0.003854
   45  H  5  Z    0.007731  -0.003853
   46  H  6  S    0.684244  -0.133288
   47  H  6  S    0.160149  -0.086283
   48  H  6  X    0.010170  -0.016435
   49  H  6  Y   -0.010352   0.019633
   50  H  6  Z   -0.010346   0.019633


                    PAIR   2

                      1          2

    1  C  1  S   -0.079049   0.182740
    2  C  1  S    0.047108  -0.488484
    3  C  1  S    0.155833  -0.060420
    4  C  1  X    0.054629  -0.292383
    5  C  1  Y   -0.066447   0.351439
    6  C  1  Z    0.066443  -0.351468
    7  C  1  X    0.044705  -0.078868
    8  C  1  Y   -0.069587   0.102943
    9  C  1  Z    0.069583  -0.102957
   10  C  1 XX    0.001182   0.010079
   11  C  1 YY   -0.000591  -0.005036
   12  C  1 ZZ   -0.000590  -0.005043
   13  C  1 XY    0.010601   0.030029
   14  C  1 XZ   -0.010601  -0.030032
   15  C  1 YZ    0.014673   0.032771
   16  H  2  S    0.684246  -0.133287
   17  H  2  S    0.160149  -0.086283
   18  H  2  X   -0.010175   0.016428
   19  H  2  Y    0.010347  -0.019633
   20  H  2  Z   -0.010346   0.019638
   21  H  3  S   -0.022537   0.039366
   22  H  3  S   -0.047925   0.018442
   23  H  3  X    0.000261  -0.004739
   24  H  3  Y   -0.007731   0.003855
   25  H  3  Z    0.007730  -0.003856
   26  C  4  S    0.012080  -0.014839
   27  C  4  S   -0.034521   0.035869
   28  C  4  S   -0.043416   0.010458
   29  C  4  X   -0.027367   0.038932
   30  C  4  Y   -0.037225   0.016190
   31  C  4  Z    0.037228  -0.016190
   32  C  4  X   -0.030075   0.018662
   33  C  4  Y   -0.016771  -0.013080
   34  C  4  Z    0.016771   0.013084
   35  C  4 XX    0.003069  -0.000606
   36  C  4 YY   -0.001534   0.000303
   37  C  4 ZZ   -0.001534   0.000303
   38  C  4 XY   -0.004902   0.006913
   39  C  4 XZ    0.004901  -0.006914
   40  C  4 YZ   -0.000879   0.000107
   41  H  5  S    0.002336  -0.004758
   42  H  5  S    0.003682  -0.007888
   43  H  5  X   -0.001081   0.000328
   44  H  5  Y   -0.000626  -0.000427
   45  H  5  Z    0.000714  -0.000278
   46  H  6  S    0.002336  -0.004757
   47  H  6  S    0.003686  -0.007890
   48  H  6  X   -0.001080   0.000328
   49  H  6  Y   -0.000714   0.000278
   50  H  6  Z    0.000626   0.000427


                    PAIR   3

                      1          2

    1  C  1  S    0.000019   0.000004
    2  C  1  S   -0.000043  -0.000010
    3  C  1  S   -0.000058  -0.000013
    4  C  1  X   -0.000019   0.000013
    5  C  1  Y    0.431378   0.025087
    6  C  1  Z    0.431392  -0.025079
    7  C  1  X   -0.000023   0.000010
    8  C  1  Y    0.385408   0.048307
    9  C  1  Z    0.385426  -0.048301
   10  C  1 XX    0.000001  -0.000001
   11  C  1 YY   -0.000001   0.000001
   12  C  1 ZZ   -0.000000  -0.000000
   13  C  1 XY   -0.000659   0.020781
   14  C  1 XZ   -0.000660  -0.020781
   15  C  1 YZ   -0.000003   0.000001
   16  H  2  S    0.000015  -0.097941
   17  H  2  S    0.000012  -0.078340
   18  H  2  X   -0.000001   0.000175
   19  H  2  Y    0.012669  -0.000791
   20  H  2  Z    0.012669   0.000792
   21  H  3  S    0.000017   0.097945
   22  H  3  S    0.000023   0.078335
   23  H  3  X   -0.000001  -0.000175
   24  H  3  Y    0.012669  -0.000792
   25  H  3  Z    0.012670   0.000791
   26  C  4  S   -0.000004  -0.000019
   27  C  4  S    0.000009   0.000040
   28  C  4  S    0.000004   0.000079
   29  C  4  X    0.000010  -0.000020
   30  C  4  Y   -0.025086  -0.431380
   31  C  4  Z   -0.025080   0.431395
   32  C  4  X   -0.000003  -0.000010
   33  C  4  Y   -0.048310  -0.385404
   34  C  4  Z   -0.048301   0.385424
   35  C  4 XX    0.000002  -0.000003
   36  C  4 YY   -0.000002   0.000002
   37  C  4 ZZ   -0.000000   0.000002
   38  C  4 XY    0.020780  -0.000659
   39  C  4 XZ    0.020781   0.000661
   40  C  4 YZ   -0.000001  -0.000004
   41  H  5  S    0.097944  -0.000012
   42  H  5  S    0.078331  -0.000005
   43  H  5  X    0.000175  -0.000001
   44  H  5  Y    0.000792  -0.012670
   45  H  5  Z    0.000792   0.012668
   46  H  6  S   -0.097945  -0.000017
   47  H  6  S   -0.078338  -0.000025
   48  H  6  X   -0.000175  -0.000000
   49  H  6  Y    0.000791  -0.012668
   50  H  6  Z    0.000791   0.012670


                    PAIR   4

                      1          2

    1  C  1  S   -0.085661   0.187356
    2  C  1  S    0.092953  -0.498274
    3  C  1  S    0.134702  -0.017607
    4  C  1  X   -0.122069   0.562741
    5  C  1  Y   -0.000003   0.000016
    6  C  1  Z    0.000001  -0.000024
    7  C  1  X   -0.108786   0.138702
    8  C  1  Y   -0.000002   0.000005
    9  C  1  Z    0.000001  -0.000005
   10  C  1 XX   -0.013831  -0.061540
   11  C  1 YY    0.006915   0.030771
   12  C  1 ZZ    0.006916   0.030770
   13  C  1 XY   -0.000005  -0.000004
   14  C  1 XZ    0.000006   0.000010
   15  C  1 YZ    0.001766  -0.003636
   16  H  2  S   -0.031053   0.039176
   17  H  2  S   -0.040796   0.015674
   18  H  2  X   -0.006047   0.005235
   19  H  2  Y    0.002956   0.002356
   20  H  2  Z   -0.002957  -0.002357
   21  H  3  S   -0.031076   0.039186
   22  H  3  S   -0.040810   0.015661
   23  H  3  X   -0.006046   0.005233
   24  H  3  Y   -0.002958  -0.002360
   25  H  3  Z    0.002958   0.002360
   26  C  4  S   -0.187352   0.085659
   27  C  4  S    0.498266  -0.092954
   28  C  4  S    0.017603  -0.134697
   29  C  4  X    0.562744  -0.122069
   30  C  4  Y   -0.000008   0.000002
   31  C  4  Z   -0.000012  -0.000001
   32  C  4  X    0.138704  -0.108785
   33  C  4  Y   -0.000002  -0.000000
   34  C  4  Z   -0.000001  -0.000001
   35  C  4 XX    0.061542   0.013832
   36  C  4 YY   -0.030771  -0.006916
   37  C  4 ZZ   -0.030771  -0.006916
   38  C  4 XY    0.000001  -0.000001
   39  C  4 XZ   -0.000008  -0.000004
   40  C  4 YZ   -0.003636   0.001766
   41  H  5  S   -0.039179   0.031062
   42  H  5  S   -0.015670   0.040801
   43  H  5  X    0.005234  -0.006046
   44  H  5  Y   -0.002357  -0.002957
   45  H  5  Z   -0.002358  -0.002957
   46  H  6  S   -0.039184   0.031073
   47  H  6  S   -0.015664   0.040809
   48  H  6  X    0.005233  -0.006046
   49  H  6  Y    0.002359   0.002958
   50  H  6  Z    0.002359   0.002958


                    PAIR   5

                      1          2

    1  C  1  S    0.012077  -0.014838
    2  C  1  S   -0.034514   0.035867
    3  C  1  S   -0.043414   0.010466
    4  C  1  X    0.027358  -0.038930
    5  C  1  Y   -0.037225   0.016189
    6  C  1  Z   -0.037230   0.016191
    7  C  1  X    0.030071  -0.018666
    8  C  1  Y   -0.016771  -0.013082
    9  C  1  Z   -0.016772  -0.013082
   10  C  1 XX    0.003069  -0.000607
   11  C  1 YY   -0.001534   0.000303
   12  C  1 ZZ   -0.001534   0.000304
   13  C  1 XY    0.004901  -0.006913
   14  C  1 XZ    0.004902  -0.006914
   15  C  1 YZ    0.000880  -0.000107
   16  H  2  S    0.002334  -0.004756
   17  H  2  S    0.003682  -0.007890
   18  H  2  X    0.001081  -0.000328
   19  H  2  Y   -0.000626  -0.000427
   20  H  2  Z   -0.000714   0.000278
   21  H  3  S    0.002335  -0.004756
   22  H  3  S    0.003683  -0.007887
   23  H  3  X    0.001080  -0.000328
   24  H  3  Y   -0.000714   0.000278
   25  H  3  Z   -0.000626  -0.000427
   26  C  4  S   -0.079050   0.182742
   27  C  4  S    0.047107  -0.488488
   28  C  4  S    0.155827  -0.060427
   29  C  4  X   -0.054626   0.292364
   30  C  4  Y   -0.066446   0.351440
   31  C  4  Z   -0.066446   0.351477
   32  C  4  X   -0.044698   0.078860
   33  C  4  Y   -0.069584   0.102947
   34  C  4  Z   -0.069588   0.102958
   35  C  4 XX    0.001186   0.010084
   36  C  4 YY   -0.000592  -0.005037
   37  C  4 ZZ   -0.000593  -0.005047
   38  C  4 XY   -0.010598  -0.030026
   39  C  4 XZ   -0.010600  -0.030031
   40  C  4 YZ   -0.014675  -0.032773
   41  H  5  S    0.684246  -0.133287
   42  H  5  S    0.160154  -0.086283
   43  H  5  X    0.010173  -0.016430
   44  H  5  Y    0.010348  -0.019633
   45  H  5  Z    0.010347  -0.019637
   46  H  6  S   -0.022544   0.039371
   47  H  6  S   -0.047925   0.018438
   48  H  6  X   -0.000260   0.004739
   49  H  6  Y   -0.007731   0.003854
   50  H  6  Z   -0.007731   0.003855


                    PAIR   6

                      1          2

    1  C  1  S   -0.079050   0.182741
    2  C  1  S    0.047111  -0.488486
    3  C  1  S    0.155825  -0.060431
    4  C  1  X    0.054620  -0.292333
    5  C  1  Y    0.066450  -0.351499
    6  C  1  Z   -0.066446   0.351444
    7  C  1  X    0.044690  -0.078839
    8  C  1  Y    0.069592  -0.102966
    9  C  1  Z   -0.069584   0.102946
   10  C  1 XX    0.001193   0.010090
   11  C  1 YY   -0.000598  -0.005051
   12  C  1 ZZ   -0.000595  -0.005039
   13  C  1 XY   -0.010597  -0.030028
   14  C  1 XZ    0.010595   0.030022
   15  C  1 YZ    0.014680   0.032777
   16  H  2  S   -0.022566   0.039379
   17  H  2  S   -0.047936   0.018428
   18  H  2  X    0.000260  -0.004741
   19  H  2  Y    0.007730  -0.003854
   20  H  2  Z   -0.007731   0.003853
   21  H  3  S    0.684244  -0.133290
   22  H  3  S    0.160156  -0.086281
   23  H  3  X   -0.010169   0.016436
   24  H  3  Y   -0.010349   0.019634
   25  H  3  Z    0.010346  -0.019632
   26  C  4  S    0.012075  -0.014838
   27  C  4  S   -0.034510   0.035869
   28  C  4  S   -0.043407   0.010471
   29  C  4  X   -0.027343   0.038925
   30  C  4  Y    0.037231  -0.016189
   31  C  4  Z   -0.037227   0.016189
   32  C  4  X   -0.030062   0.018670
   33  C  4  Y    0.016774   0.013082
   34  C  4  Z   -0.016772  -0.013082
   35  C  4 XX    0.003068  -0.000607
   36  C  4 YY   -0.001534   0.000304
   37  C  4 ZZ   -0.001534   0.000304
   38  C  4 XY    0.004902  -0.006914
   39  C  4 XZ   -0.004901   0.006913
   40  C  4 YZ   -0.000879   0.000108
   41  H  5  S    0.002332  -0.004754
   42  H  5  S    0.003681  -0.007889
   43  H  5  X   -0.001080   0.000328
   44  H  5  Y    0.000714  -0.000278
   45  H  5  Z   -0.000626  -0.000427
   46  H  6  S    0.002331  -0.004751
   47  H  6  S    0.003676  -0.007881
   48  H  6  X   -0.001080   0.000328
   49  H  6  Y    0.000626   0.000427
   50  H  6  Z   -0.000714   0.000278

 ... END OF ROHF-GVB SCF CALCULATION ...
 STEP CPU TIME =     6.86 TOTAL CPU TIME =          6.9 (      0.1 MIN)
 TOTAL WALL CLOCK TIME=          7.5 SECONDS, CPU UTILIZATION IS    92.62%

     ----------------------------------------------------------------
     PROPERTY VALUES FOR THE GVB   SELF-CONSISTENT FIELD WAVEFUNCTION
     ----------------------------------------------------------------

          -----------------
          ENERGY COMPONENTS
          -----------------

         WAVEFUNCTION NORMALIZATION =       1.0000000000

                ONE ELECTRON ENERGY =    -166.5548946379
                TWO ELECTRON ENERGY =      56.8868102780
           NUCLEAR REPULSION ENERGY =      31.6688684670
                                      ------------------
                       TOTAL ENERGY =     -77.9992158929

 ELECTRON-ELECTRON POTENTIAL ENERGY =      56.8868102780
  NUCLEUS-ELECTRON POTENTIAL ENERGY =    -244.1167349668
   NUCLEUS-NUCLEUS POTENTIAL ENERGY =      31.6688684670
                                      ------------------
             TOTAL POTENTIAL ENERGY =    -155.5610562218
               TOTAL KINETIC ENERGY =      77.5618403289
                 VIRIAL RATIO (V/T) =       2.0056390560

          ---------------------------------------
          MULLIKEN AND LOWDIN POPULATION ANALYSES
          ---------------------------------------

     ATOMIC MULLIKEN POPULATION IN EACH MOLECULAR ORBITAL

                      1          2          3          4          5

                  2.000000   2.000000   1.985420   0.014580   1.985419

    1             1.004904   0.994758  -0.007409   0.000096   1.063465
    2             0.000263  -0.000094   0.000579   0.000001   0.941571
    3             0.000263  -0.000094   0.000578   0.000001  -0.013365
    4             0.994043   1.005619   1.063465   0.007635  -0.007410
    5             0.000263  -0.000094  -0.013362  -0.000009   0.000579
    6             0.000263  -0.000094   0.941570   0.006856   0.000579

                      6          7          8          9         10

                  0.014581   1.000052   0.999948   1.989313   0.010687

    1             0.007635   0.482787   0.482736   1.017011   0.005322
    2             0.006857   0.008618   0.008620  -0.011178   0.000011
    3            -0.000009   0.008621   0.008618  -0.011174   0.000011
    4             0.000096   0.482787   0.482736   1.017005   0.005322
    5             0.000001   0.008618   0.008619  -0.011176   0.000011
    6             0.000001   0.008621   0.008618  -0.011175   0.000011

                     11         12         13         14

                  1.985419   0.014581   1.985420   0.014580

    1            -0.007412   0.000096   1.063459   0.007635
    2             0.000579   0.000001  -0.013361  -0.000009
    3             0.000579   0.000001   0.941575   0.006856
    4             1.063462   0.007635  -0.007411   0.000096
    5             0.941575   0.006857   0.000579   0.000001
    6            -0.013363  -0.000009   0.000578   0.000001

               ----- POPULATIONS IN EACH AO -----
                             MULLIKEN      LOWDIN
              1  C  1  S      1.98997     1.94529
              2  C  1  S      0.88592     0.42077
              3  C  1  S      0.29322     0.26708
              4  C  1  X      0.71630     0.63989
              5  C  1  Y      0.63154     0.57981
              6  C  1  Z      0.63153     0.57980
              7  C  1  X      0.21591     0.35073
              8  C  1  Y      0.36114     0.41988
              9  C  1  Z      0.36115     0.41988
             10  C  1 XX      0.00223     0.18343
             11  C  1 YY      0.00669     0.13686
             12  C  1 ZZ      0.00812     0.13686
             13  C  1 XY      0.00812     0.01923
             14  C  1 XZ      0.00325     0.01923
             15  C  1 YZ      0.00000     0.00882
             16  H  2  S      0.75240     0.62187
             17  H  2  S      0.16635     0.25189
             18  H  2  X      0.00606     0.01575
             19  H  2  Y      0.00882     0.02335
             20  H  2  Z      0.00882     0.02336
             21  H  3  S      0.75240     0.62187
             22  H  3  S      0.16636     0.25189
             23  H  3  X      0.00606     0.01575
             24  H  3  Y      0.00882     0.02336
             25  H  3  Z      0.00882     0.02335
             26  C  4  S      1.98997     1.94529
             27  C  4  S      0.88592     0.42077
             28  C  4  S      0.29321     0.26708
             29  C  4  X      0.71630     0.63988
             30  C  4  Y      0.63154     0.57981
             31  C  4  Z      0.63153     0.57980
             32  C  4  X      0.21591     0.35073
             33  C  4  Y      0.36114     0.41988
             34  C  4  Z      0.36115     0.41988
             35  C  4 XX      0.00223     0.18343
             36  C  4 YY      0.00669     0.13686
             37  C  4 ZZ      0.00812     0.13686
             38  C  4 XY      0.00812     0.01923
             39  C  4 XZ      0.00325     0.01923
             40  C  4 YZ      0.00000     0.00882
             41  H  5  S      0.75240     0.62187
             42  H  5  S      0.16636     0.25189
             43  H  5  X      0.00606     0.01575
             44  H  5  Y      0.00882     0.02335
             45  H  5  Z      0.00882     0.02336
             46  H  6  S      0.75240     0.62187
             47  H  6  S      0.16635     0.25189
             48  H  6  X      0.00606     0.01575
             49  H  6  Y      0.00882     0.02336
             50  H  6  Z      0.00882     0.02335

          ----- MULLIKEN ATOMIC OVERLAP POPULATIONS -----
          (OFF-DIAGONAL ELEMENTS NEED TO BE MULTIPLIED BY 2)

             1           2           3           4           5

    1    5.0614512
    2    0.4080063   0.6016963
    3    0.4080038  -0.0305434   0.6016996
    4    0.3154870  -0.0389358  -0.0389301   5.0614477
    5   -0.0389349   0.0011168   0.0011169   0.4080060   0.6016994
    6   -0.0389308   0.0011171   0.0011147   0.4080051  -0.0305429

             6

    6    0.6016947

          TOTAL MULLIKEN AND LOWDIN ATOMIC POPULATIONS
       ATOM         MULL.POP.    CHARGE          LOW.POP.     CHARGE
    1 C             6.115083   -0.115083         6.127553   -0.127553
    2 H             0.942457    0.057543         0.936223    0.063777
    3 H             0.942461    0.057539         0.936225    0.063775
    4 C             6.115080   -0.115080         6.127552   -0.127552
    5 H             0.942461    0.057539         0.936224    0.063776
    6 H             0.942458    0.057542         0.936223    0.063777

          MULLIKEN SPHERICAL HARMONIC POPULATIONS
       ATOM           S       P       D      F      G      H      I    TOTAL
    1 C             3.17    2.92    0.03   0.00   0.00   0.00   0.00    6.12
    2 H             0.92    0.02    0.00   0.00   0.00   0.00   0.00    0.94
    3 H             0.92    0.02    0.00   0.00   0.00   0.00   0.00    0.94
    4 C             3.17    2.92    0.03   0.00   0.00   0.00   0.00    6.12
    5 H             0.92    0.02    0.00   0.00   0.00   0.00   0.00    0.94
    6 H             0.92    0.02    0.00   0.00   0.00   0.00   0.00    0.94

          ---------------------
          ELECTROSTATIC MOMENTS
          ---------------------

 POINT   1           X           Y           Z (BOHR)    CHARGE
                 0.000000   -0.000010   -0.000000        0.00 (A.U.)
         DX          DY          DZ         /D/  (DEBYE)
    -0.000010    0.000149    0.000014    0.000150
 ...... END OF PROPERTY EVALUATION ......
 STEP CPU TIME =     0.01 TOTAL CPU TIME =          6.9 (      0.1 MIN)
 TOTAL WALL CLOCK TIME=          7.5 SECONDS, CPU UTILIZATION IS    92.63%
                580000  WORDS OF DYNAMIC MEMORY USED
 EXECUTION OF GAMESS TERMINATED NORMALLY Wed Jan 20 01:39:27 2021
 DDI: 263640 bytes (0.3 MB / 0 MWords) used by master data server.

 ----------------------------------------
 CPU timing information for all processes
 ========================================
 0: 6.737 + 0.180 = 6.918
 ----------------------------------------
 ddikick.x: exited gracefully.
unset echo
----- accounting info -----
Files used on the master node xn01 were:
-rw-rw-r-- 1 srwang srwang   42736 Jan 20 01:39 /scratch/scr/srwang/gamess/eth_uhf_uno_asrot2gvb6.dat
-rw-rw-r-- 1 srwang srwang   42242 Jan 20 01:39 /scratch/scr/srwang/gamess/eth_uhf_uno_asrot2gvb6.F05
-rw-rw-r-- 1 srwang srwang 4515360 Jan 20 01:39 /scratch/scr/srwang/gamess/eth_uhf_uno_asrot2gvb6.F10
-rw-rw-r-- 1 srwang srwang  312984 Jan 20 01:39 /scratch/scr/srwang/gamess/eth_uhf_uno_asrot2gvb6.F15
-rw-rw-r-- 1 srwang srwang  387904 Jan 20 01:39 /scratch/scr/srwang/gamess/eth_uhf_uno_asrot2gvb6.F23
ls: No match.
ls: No match.
ls: No match.
Wed Jan 20 01:39:30 CST 2021
0.209u 0.126s 0:10.82 2.9%	0+0k 0+88io 0pf+0w
