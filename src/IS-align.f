**************************************************************************
*     iAlign is a program for assessing protein-protein interface        *
*     similarity. The source code evolved from TM-align.                 *
*                                                                        *
*     Mu Gao  <mu.gao.@gatech.edu>                                       *
**************************************************************************
      
      program isalign
      include 'pars.h'

      COMMON/BACKBONE/XA(3,maxr,0:1)  !coordinates of ca of two structures for comparison
      integer invmap0(maxr)     !final alignment from second sequence to first
      integer invmap0_r(maxr)   !final alignment from first  sequence to second
      common/length/nseq1,nseq2 !length of structures (number of residues)
      common/d0/d0,anseq
      common/d0min/d0_min
      common/d00/d00,d002
      common/d8/d8

      character*500 fnam,pdb(100),outname,pdb1,pdb2,file1,file2
      character*500 contfile1,contfile2    !contact list files
      character*3   ss1(maxr),ss2(maxr)    !three letter AA name 
      character     seq1(maxr),seq2(maxr)  !one letter AA name
      character     aseq1(maxr),aseq2(maxr),aseq3(maxr)  !aligned sequences
      character*80  matchlst(maxr)    !list of matched residues
      character*2   measure   ! 'TM' -> TM-score, 'IS' -> IS-score

      integer   m1(maxr),m2(maxr)
      integer   mf1(maxr),mf2(maxr)
      real      xtm1(maxr),ytm1(maxr),ztm1(maxr)
      real      xtm2(maxr),ytm2(maxr),ztm2(maxr)
      real      xtmf1(maxr),ytmf1(maxr),ztmf1(maxr)
      real      xtmf2(maxr),ytmf2(maxr),ztmf2(maxr)

ccc   RMSD:
      double precision r_1(3,maxr),r_2(3,maxr),r_3(3,maxr),w(maxr)
      double precision u(3,3),t(3),rms,drms !armsd is real
      double precision u_final(3,3),t_final(3)
      data w /maxr*1.0/

      common/id2chain/id_chain1(maxr),id_chain2(maxr)  !original index to chain ID
      integer cont1(maxc,maxr),cont2(maxc,maxr)     !contact list
      common/contlst/cont1,cont2
      common/contnum/ncont1(maxr),ncont2(maxr)      !number of contacts
      common/contcf/d_col,d_col2

      common/chains/ichainterm1(0:maxk),ichainterm2(0:maxk) !indexes of C-termini of chains
      common/sec/isec(maxr),jsec(maxr)   !secondary structure

      real   ca_ca_cf   ! distance cutoff for defining Calpha-Calpha contact 

      character*5 respdbid1(maxr),respdbid2(maxr)
      character   chname1(maxk),  chname2(maxk)
      character   pdbchnm1(maxk), pdbchnm2(maxk)
      integer  nchain1,nchain2,pdbnchain1,pdbnchain2
      integer  lstr1,lstr2,nseq


ccc   Timing:
      real*8 total_time
      integer ( kind = 4 ) clock_start
      integer ( kind = 4 ) clock_stop
      integer ( kind = 4 ) clock_max
      integer ( kind = 4 ) clock_rate



ccc   Options:
      integer  score_flag  ! 0 - TM-score, 1 - IS-score
      integer  order_flag  ! 0 - sequential order independent, 1 - dependent
      integer  pv_flag     ! 0 - use score based on short sequence, 1 - regular
      integer  quick_flag  ! 0 - regular score caculations, 1 - fast score calculations (default), 2 - two rounds
      integer  verbo_flag  ! 2 - print detailed match list, default 1 in sequential, 2 in non-sequential mode
      integer  cont_flag   ! 0 - read contact from file, 1 - calculate interface Ca-Ca contacts
      integer  seqd0_flag   ! 0 - use non-sequential d0, 1 - use sequential d0 (for testing purpose)

      integer  tarlen, temlen  ! target and template length for p-value calculations
      integer  best_iter, best_init, best_nal, best_ncol, lseq,mode
      real     max_TM, max_TM8, min_RM8, diff_TM8, diff_RM8
      real     best_rmsd, best_sid, scalelen,d0_final,f0
      real*8   pvalue, best_pvalue, zs, best_zs
      real     eps


      call system_clock ( clock_start, clock_rate, clock_max )


c================================================================================
      call getarg(1,fnam)
      if(fnam.eq.' '.or.fnam.eq.'?'.or.fnam.eq.'-h')then
         call printHelp()
         goto 9999
      endif
 
******* options ----------->
      m_out=-1            !decided output
      m_fix=-1            !fixed length-scale only for output
      m_ave=-1            !using average length
      m_d0_min=-1         !diminum d0 for search
      m_d0=-1             !given d0 for both search and output
      score_flag=1        !default: 1 - IS-score, 0 - TM-score
      measure='IS'
      order_flag=1        !default: 1 - order dependent, 0 - order independent
      pv_flag=0           !default: 0 - regular, 3 - no pv calculation
      quick_flag=1
      verbo_flag=1
      seqd0_flag=0
      cont_flag=0
      pvalue=-1
      zscore=0
      eps=1.e-6
      narg=iargc()
      i=0
      j=0
 115  continue
      i=i+1
      call getarg(i,fnam)
      if(fnam.eq.'-o')then
         m_out=1
         i=i+1
         call getarg(i,outname)
      elseif(fnam.eq.'-L')then  !change both L_all and d0
         m_fix=1
         i=i+1
         pv_flag=3
         call getarg(i,fnam)
         read(fnam,*)L_fix
      elseif(fnam.eq.'-dmin')then
         m_d0_min=1
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)d0_min_input
         pv_flag=3
      elseif(fnam.eq.'-d0')then
         m_d0=1
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)d0_fix
      elseif(fnam.eq.'-a')then !this will change TM-score but not the alignment
         m_ave=1
         pv_flag=3
      elseif(fnam.eq.'-b')then
         m_ave=2
      elseif(fnam.eq.'-c')then
         m_ave=3
      elseif(fnam.eq.'-t')then ! use regular TM-score measure
         score_flag=0
         measure='TM'
      elseif(fnam.eq.'-s')then ! sequential-order-independent
         order_flag=0
         if(verbo_flag.ne.0)verbo_flag=2
      elseif(fnam.eq.'-ca')then ! calculate contacts
         cont_flag=1
      elseif(fnam.eq.'-seqd0')then ! use sequential d0 for testing purpose 
         seqd0_flag=1
      else if(fnam.eq.'-v')then ! print detailed match list
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)verbo_flag
         if(verbo_flag.lt.0.or.verbo_flag.gt.2)then
            call printHelp()
            goto 9999
         endif
      elseif(fnam.eq.'-q')then ! quick mode
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)quick_flag
         if(quick_flag.lt.0.or.quick_flag.gt.2)then
            call printHelp()
            goto 9999
         endif
      else
         j=j+1
         pdb(j)=fnam
      endif
      if(i.lt.narg)goto 115

      if(order_flag.eq.0.and.verbo_flag.eq.1) verbo_flag=2
      if(order_flag.eq.0.and.score_flag.eq.0) pv_flag=3  ! TM-score is not informative in non-sequential mode


ccccc read PDB files
      call readPDBFile(pdb(1),xa,0,nchain1,ichainterm1,pdbchnm1,
     &     seq1,ss1,respdbid1,isec,lstr1,cont_flag,cont1,ncont1)
      temlen=lstr1
      nseq1=lstr1
      do i=1,nchain1
         chname1(i)=pdbchnm1(i)
      enddo

      call readPDBFile(pdb(2),xa,1,nchain2,ichainterm2,pdbchnm2,
     &     seq2,ss2,respdbid2,jsec,lstr2,cont_flag,cont2,ncont2)
      tarlen=lstr2
      nseq2=lstr2
      do i=1,nchain2
         chname2(i)=pdbchnm2(i)
      enddo


ccccc obtain interface contacts
      if(cont_flag.eq.0)then
         contfile1=pdb(3)
         contfile2=pdb(4)
         call readContFile(contfile1,ncont1,cont1)  !read contact lists pre-calculated
         call readContFile(contfile2,ncont2,cont2)
      endif


      !call timestamp ( )
cccccc Error message in case there is only 1 chain 

      if (nseq1.le.5.or.nseq2.le.5)then
         write(*,*)
         write(*,*)'Error: the minimum length of a structure'//
     &        ' is six residues.'
         write(*,*)
         goto 9999
      endif

      
ccccc---------------- set parameters here -------------------------cccccc

*!!!  Scale of TM-score in search is based on the smaller protein --------->
      d0_min=0.5
      if(m_d0_min.eq.1)then
         d0_min=d0_min_input    !for search
      endif
      anseq_min=min(nseq1,nseq2)
      nseq=max(nseq1,nseq2)
      anseq=anseq_min           !length for defining TMscore in search
      if(anseq.gt.15)then
         d0=1.24*(anseq-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
         d0=d0_min
      endif

      if(order_flag.eq.0.and.seqd0_flag.eq.0)then
         if(anseq.gt.15)then
            d0=0.7*(anseq-15)**(1./3.) - 0.1
         else
            d0=0.6
         endif
      endif

      if(d0.lt.d0_min)d0=d0_min
      if(m_d0.eq.1)d0=d0_fix

      d8=1.5*anseq_min**0.3+3.5 !remove pairs with dis>d8 during search & final
      if(order_flag.eq.0.and.seqd0_flag.eq.0)then
         if(d8.gt.8)d8=8
      endif

      !scale length for calculating final TMscore
      scalelen=nseq2               
      if(m_ave.eq.1)scalelen=(nseq1+nseq2)/2.0 !<L>
      if(m_ave.eq.2)scalelen=min(nseq1,nseq2)
      if(m_ave.eq.3)scalelen=max(nseq1,nseq2)
      if(scalelen.lt.anseq_min)scalelen=anseq_min
      if(m_fix.eq.1)scalelen=L_fix !input length
      if(scalelen.gt.15)then
         d0_final=1.24*(scalelen-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
         d0_final=d0_min
      endif

      if(order_flag.eq.0.and.seqd0_flag.eq.0)then
         if(scalelen.gt.15)then
            !d0_final=0.55*(scalelen-20)**(1./3.)+0.3
            d0_final=0.7*(scalelen-15)**(1./3.) - 0.1
         else
            d0_final=0.6
         endif
      endif

      if(d0_final.lt.d0_min)d0_final=d0_min
      if(m_d0.eq.1)d0_final=d0_fix

      d00=d0                    !for quickly calculate TM-score in searching
      if(d00.gt.8)d00=8
      if(d00.lt.4.5)d00=4.5
      d002=d00**2

      d_col=d8          !distance cutoff for aligned residues with contact overlap
      d_col2=d_col**2   !prevent costs of sqrt

      d_output=5.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      num_cont1=0  !number of interfacial contacts of protein 1
      num_cont2=0  !number of interfacial contacts of protein 2
      do i=1,nseq1
         if(i.le.ichainterm1(1))then
            id_chain1(i)=1
            num_cont1=num_cont1 + ncont1(i)  ! only need count once
         else
            id_chain1(i)=2
         endif
      enddo
      do j=1,nseq2
         if(j.le.ichainterm2(1))then
            id_chain2(j)=1
            num_cont2=num_cont2 + ncont2(j)
         else
            id_chain2(j)=2
         endif
      enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(num_cont1.lt.5.or.num_cont2.lt.5)then
         write(*,'(2(A,I3))')'num_cont1 =',num_cont1,', num_con2=',num_cont2
         write(*,'(A)')'Error: the number of contacts is too small.'//
     &        ' Please check your input files.'
         write(*,*)
         goto 9999
      endif

      if(min(nseq1,nseq2).lt.15)then
         pv_flag=3              !no p-value calculation for very short interface
         write(*,'(A)')'Warning: a minimum of 15 residues is required '//
     &        'for p-value calculations. P-value not caculated!'
      endif


      max_TM=-1
      max_TM8=-1
      min_RM8=100
      best_iter=1
      best_rmsd=0
      best_nal=0
      best_init=0
      best_ncol=0
      best_pvalue=1
      best_zs=0
ccccc----------------------------------------------------------------------------
ccccc Start the main loop of swapping chains. This solves the issue 
ccccc of score-dependence on the order of chains. no need to swap structures
ccccc----------------------------------------------------------------------------
      do 8000 iter=1,4
         if(max_TM8.gt.0.99.and.iter.gt.1) goto 8001  !jump out of the loop
         if(max_TM8.gt.0.9 .and.iter.gt.2) goto 8001  !jump out of the loop
         if(quick_flag.eq.2.and.iter.gt.2) goto 8001  !jump out of the loop

         !write(*,'(A,I2)') '**** Round: ',iter
         if(iter.eq.2)then
            call swap_chains(nchain2,ichainterm2,chname2,seq2,xa,1,
     &           ss2,respdbid2,jsec,id_chain2,ncont2,cont2)

         else if(iter.eq.3)then
            call swap_chains(nchain1,ichainterm1,chname1,seq1,xa,0,
     &           ss1,respdbid1,isec,id_chain1,ncont1,cont1)

         else if(iter.eq.4)then
            call swap_chains(nchain2,ichainterm2,chname2,seq2,xa,1,
     &           ss2,respdbid2,jsec,id_chain2,ncont2,cont2)

         else if(iter.eq.5)then
            call swap_structs(xa,nchain1,ichainterm1,chname1,seq1,0,
     &           ss1,respdbid1,isec,id_chain1,ncont1,cont1,nseq1,
     &           nchain2,ichainterm2,chname2,seq2,1,
     &           ss2,respdbid2,jsec,id_chain2,ncont2,cont2,nseq2)

         else if(iter.eq.6)then
            call swap_chains(nchain2,ichainterm2,chname2,seq2,xa,1,
     &           ss2,respdbid2,jsec,id_chain2,ncont2,cont2)

         else if(iter.eq.7)then
            call swap_chains(nchain1,ichainterm1,chname1,seq1,xa,0,
     &           ss1,respdbid1,isec,id_chain1,ncont1,cont1)

         else if(iter.eq.8)then
            call swap_chains(nchain2,ichainterm2,chname2,seq2,xa,1,
     &           ss2,respdbid2,jsec,id_chain2,ncont2,cont2)
         endif


***** find alignment from seq2 to seq1 **************************

         call super_align(invmap0,nseq,nseq1,nseq2,score_flag,
     &        order_flag,quick_flag,init,mode)

*********************************************

        !resuperpose to find residues of dis<d8 ------------------------>
         n_al=0
         if(mode.eq.0)then
            do j=1,nseq2
               if(invmap0(j).gt.0)then
                  i=invmap0(j)
                  n_al=n_al+1
                  xtm1(n_al)=xa(1,i,0)
                  ytm1(n_al)=xa(2,i,0)
                  ztm1(n_al)=xa(3,i,0)
                  xtm2(n_al)=xa(1,j,1)
                  ytm2(n_al)=xa(2,j,1)
                  ztm2(n_al)=xa(3,j,1)
                  m1(n_al)=i    !for recording residue order
                  m2(n_al)=j
               endif
            enddo
         else
            do i=1,nseq1
               invmap0_r(i)=-1   !seq1 to seq2, reverse of invmap
            enddo

            do j=1,nseq2
               i=invmap0(j)
               if(i.gt.0) invmap0_r(i)=j
            enddo

            do i=1,nseq1
               j=invmap0_r(i)    !j aligned to i
               if(j.gt.0)then
                  n_al=n_al+1
                  xtm1(n_al)=xa(1,i,0) !for TM-score
                  ytm1(n_al)=xa(2,i,0)
                  ztm1(n_al)=xa(3,i,0)
                  xtm2(n_al)=xa(1,j,1)
                  ytm2(n_al)=xa(2,j,1)
                  ztm2(n_al)=xa(3,j,1)

                  m1(n_al)=i !aligned position to original index
                  m2(n_al)=j
               endif
            enddo
         endif
          
         d0_input=d0
         d0_search=d0
         if (d0_search.gt.8.0)d0_search=8.0
         if (d0_search.lt.3.0)d0_search=4.5
         call TMsearch(d0_input,d0_search,n_al,xtm1,ytm1,ztm1,m1,n_al,
     &        xtm2,ytm2,ztm2,m2,TM,Rcomm,Lcomm,2,score_flag,ncol) !calcualte TMscore with dis<d8 only
         TM8=TM*n_al/anseq
         !print *,'TMscore8 ',TM8  !TM-score8
      
         !remove dis>d8 in normal TM-score calculation for final report----->
         j=0
         n_eq=0
         do i=1,n_al
            dis2=sqrt((xtm1(i)-xtm2(i))**2+(ytm1(i)-ytm2(i))**2+
     &           (ztm1(i)-ztm2(i))**2)
            if(dis2.le.d8)then
               j=j+1
               xtm1(j)=xtm1(i)
               ytm1(j)=ytm1(i)
               ztm1(j)=ztm1(i)
               xtm2(j)=xtm2(i)
               ytm2(j)=ytm2(i)
               ztm2(j)=ztm2(i)
               m1(j)=m1(i)
               m2(j)=m2(i)
               if(ss1(m1(i)).eq.ss2(m2(i)))then
                  n_eq=n_eq+1
               endif
               !print *,j,' m1=',m1(j),' m2=',m2(j)
            endif
         enddo
         
         n8_al=j
         seq_id=float(n_eq)/float(n8_al)   !sequence identity over aligned region


         !for calculating p-value of scores, always use short sequence length to derive d0_search
         !recalculate the score is necessary since aligned residues are slightly different by excluding d > d8
         call TMsearch(d0_input,d0_search,n8_al,xtm1,ytm1,ztm1,
     &        m1,n8_al,xtm2,ytm2,ztm2,m2,TM,Rcomm,Lcomm,3,
     &        score_flag,ncol)
         TM8=TM*n8_al/anseq
         RM8=Rcomm

         if(pv_flag.le.2)then
            if( score_flag.eq.1 )then
               f0=0.18 - 0.35*(anseq**(-0.3))
               TM8=(TM8+f0)/(1+f0)       ! rescale IS-score
            endif
            call calcPvalue(TM8,score_flag,order_flag,temlen,tarlen,pvalue,zs)
         endif

         !final score
         d0_input=d0_final
         d0_search=d0_final
         if (d0_search.gt.8.0)d0_search=8.0
         if (d0_search.lt.4.5)d0_search=4.5

         call TMsearch(d0_input,d0_search,n8_al,xtm1,ytm1,ztm1,m1,
     &        n8_al,xtm2,ytm2,ztm2,m2,TM,Rcomm,Lcomm,3,score_flag,ncol) !normal TMscore
         TM=TM*n8_al/scalelen
         if(score_flag.eq.1)then
            f0 = 0.18 - 0.35*(scalelen**(-0.3))
            TM=(TM+f0)/(1+f0)   ! rescale IS-score
         endif

!         write(*,'(A6,I2,A12,F10.7,A5,F10.7,A9,F10.7)')'Iter=',iter,
!     &        ' fTMscore =',TM, ' TM8=',TM8, ' max_TM8=',max_TM8
!         write(*,'(A5,F8.5,A9,F9.5, /)')' RM8= ',Rcomm,' min_RM8=',min_RM8

         ! save the best alignment. To avoid possible asymmetry of alignment,
         ! always use the best alignment obtained with the shorter sequence.
         diff_TM8 = TM8 - max_TM8
         diff_RM8 = min_RM8 - RM8
         if((diff_TM8.gt.eps).or.
     &        (abs(diff_TM8).lt.eps.and.diff_RM8.gt.eps))then
            do i=1,n8_al
               xtmf1(i)=xtm1(i)
               ytmf1(i)=ytm1(i)
               ztmf1(i)=ztm1(i)
               xtmf2(i)=xtm2(i)
               ytmf2(i)=ytm2(i)
               ztmf2(i)=ztm2(i)
               mf1(i)=m1(i)
               mf2(i)=m2(i)               
            enddo

            if(score_flag.eq.0)call cal_col(mf1,mf2,n8_al,ncol) !number of contact overlaps

            max_TM     = TM
            max_TM8    = TM8
            min_RM8    = RM8
            best_nal   = n8_al
            best_rmsd  = Rcomm
            best_sid   = seq_id
            best_ncol  = ncol
            best_iter  = iter
            best_init  = init
            best_pvalue= pvalue
            best_zs    = zs


            ! save the best transformation matrix
            call get_TransMatrix(mf1,xa,xtmf1,ytmf1,ztmf1,
     &           best_nal,t_final,u_final)

            ! save the best sequential alignment

            if(verbo_flag.eq.2) then
               call get_matchlst(nseq1,nseq2,mf1,mf2,xtmf1,xtmf2,
     &              ytmf1,ytmf2,ztmf1,ztmf2,respdbid1,respdbid2,
     &              seq1,seq2,id_chain1,id_chain2,ncont1,ncont2,
     &              chname1,chname2,best_nal,matchlst,best_iter)
            else if(verbo_flag.eq.1) then
               call get_alnseq(mf1,mf2,best_nal,seq1,seq2,nchain1,
     &              nchain2,ichainterm1,ichainterm2,xtmf1,ytmf1,
     &              ztmf1,xtmf2,ytmf2,ztmf2,d_output,aseq1,aseq2,
     &              aseq3,lseq)
            endif
         endif

 8000 continue                  !enddo
 8001 continue
********* Score calculations are done! *******************


********* Output results ******************************
      write(*,'(/,2A)')'***********************************************',
     &     '**********************'
      write(*,'(2A)')'                             iAlign           '
      write(*,'(2A)')'  A tool for structural comparison of ',
     &     'protein-protein interfaces.'
      write(*,'(2A)')'  Please report bugs to:  Mu Gao  ',
     &     '<mu.gao@gatech.edu>.        '
      write(*,'(2A,/)')'***********************************************',
     &     '**********************'

      call getBasename(pdb(1),file1)
      call getBasename(pdb(2),file2)

      write(*,101)1,file1,pdbchnm1(1),pdbchnm1(2),lstr1,num_cont1
 101  format('Structure ',I1,': ',A20,' Chains ',A1,1X,A1,','
     &     I4, ' AAs,',I4,' Contacts')

      write(*,101)2,file2,pdbchnm2(1),pdbchnm2(2),lstr2,num_cont2
      write(*,*)

      if(best_nal.lt.3)then
         write(*,*)
         write(*,*)'No significant alignment found!'
         write(*,*)
         goto 9999
      endif

      if(pv_flag.le.2)then
         write(*,103)measure,max_TM,best_pvalue,best_zs
         write(*,105)best_nal,best_ncol,best_rmsd,best_sid
      else  ! no p-value output
         write(*,104)measure,max_TM
         write(*,105)best_nal,best_ncol,best_rmsd,best_sid
      endif

 103  format(A2,'-score =',F8.5,', P-value =',E12.4E3,', Z-score =',F7.3)
 104  format(A2,'-score =',F8.5)
 105  format('Number of aligned residues  =',I4, /,
     &       'Number of aligned contacts  =',I4, /,
     &       'RMSD =',F6.2, ', Seq identity  =',F6.3, /)



********* print the transformation matrix ------------>
      !if(best_nal.le.3) goto 8888  !skip matrix

      if(best_iter.ge.5)then ! proteins were swapped
         call inv_TransMatrix( t_final, u_final )
      endif

      write(*,*)'----- Transformation matrix for aligning ',
     &     'struct 1 to struct 2 -----'
      write(*,*)'i          t(i)         r(i,1)         r(i,2) ',
     &     '        r(i,3)'
      do i=1,3
         write(*,204)i,t_final(i),u_final(i,1),u_final(i,2),u_final(i,3)
      enddo
      write(*,*)
 204  format(I2,f18.10,f15.10,f15.10,f15.10)

 8888 continue
      

********* print the alignment **************************
      if( verbo_flag.eq.1 ) then
         write(*,107) d_output
 107     format('Interface Alignment (":" AA pairs within ',F3.1,
     &        ' Angstrom in Ca-Ca distance)')

         if     (best_iter.eq.1.or.best_iter.eq.6)then
            write(*,108)pdbchnm1(1),pdbchnm1(2),pdbchnm2(1),pdbchnm2(2)

         else if(best_iter.eq.2.or.best_iter.eq.7)then
            write(*,108)pdbchnm1(1),pdbchnm1(2),pdbchnm2(2),pdbchnm2(1)

         else if(best_iter.eq.3.or.best_iter.eq.8)then
            write(*,108)pdbchnm1(2),pdbchnm1(1),pdbchnm2(2),pdbchnm2(1)

         else if(best_iter.eq.4.or.best_iter.eq.5)then
            write(*,108)pdbchnm1(2),pdbchnm1(1),pdbchnm2(1),pdbchnm2(2)
         endif
 108     format('Struct 1 Chains ',A1,A1,1X,' vs. ',
     &        'Struct 2 Chains ',A1,A1,
     &        ' (AAs in upper/lower cases)',/)

         if(best_iter.lt.5)then !proteins were not swapped
            write(*,10)(aseq1(i),i=1,lseq)
            write(*,10)(aseq3(i),i=1,lseq)
            write(*,10)(aseq2(i),i=1,lseq)
         else
            write(*,10)(aseq2(i),i=1,lseq)
            write(*,10)(aseq3(i),i=1,lseq)
            write(*,10)(aseq1(i),i=1,lseq)
         endif
 10      format(10000A1)

         write(*,*)
      endif
      if( verbo_flag.eq.2 ) then
         write(*,*)'-----   Aligned Interface Residues  -----'
         write(*,*)'----- Structure 1       Structure 2 -----'
         write(*,*)'Index Ch1 Resid1 AA1    Ch2 Resid2 AA2   ',
     &     ' Distance NAC NC1 NC2 Note'
         do i = 1, best_nal
            write(*,'(A70)') matchlst(i)
         enddo

         write(*,*)'":" AA pairs within 5 Ansgtrom in Ca-Ca distance'
         write(*,*)'"*" Identical AA pairs'
         write(*,*)
      endif

      write(*,109)nint(scalelen),d0_final
 109  format('Scoring parameters: normalization length =',I4,',  d0 =',F6.3)

      if(order_flag.eq.0) then
         write(*,'(A)')'Alignment search mode: non-sequential'
      else
         write(*,'(A)')'Alignment search mode: sequential'
      endif

      write(*,110)best_iter,best_init
 110  format('Best alignment search: round =',I2,
     &     ', initial =',I3)

      !call timestamp ( )

      call system_clock ( clock_stop, clock_rate, clock_max )
      total_time = (clock_stop - clock_start)/real(clock_rate, kind=8)

      write(*,'(A,F12.5,A,/)')'Running time: ', total_time, ' seconds'

 9999 END

***---------------------  Main procedure ends here ----------------***





ccc================================================================cccc
ccc===============         Subroutines     ========================cccc 
ccc================================================================cccc


***********************************************************************
***********************************************************************
*     Structure superposition
***********************************************************************
***********************************************************************
      SUBROUTINE super_align(invmap0,nseq,nseq1,nseq2,score_flag,
     &     order_flag,quick_flag,best_init,best_mode)
      implicit none
      include 'pars.h'

      common/id2chain/id_chain1(maxr),id_chain2(maxr)
      common/contlst/cont1(maxc,maxr),cont2(maxc,maxr) !contact list
      common/contnum/ncont1,ncont2       !number of contacts

      integer id_chain1, id_chain2
      integer ncont1(maxr),ncont2(maxr)             
      integer cont1,cont2     

      COMMON/BACKBONE/XA
      real    XA(3,maxr,0:1),coor1(3,maxr),coor2(3,maxr)

      common/chains/ichainterm1,ichainterm2 !indexes of C-termini of chains
      common/sec/isec,jsec   !secondary structure

      integer nseq1,nseq2,nseq
      integer len1,len2,len11,len12,len21,len22,maxl
      integer ichainterm1(0:maxk),ichainterm2(0:maxk)
      integer isec(maxr),jsec(maxr)

      common/TM/TM,TMmax
      real    TM,TMmax,BTMmax,TM_old
      real    TMcut1,TMcut2  !cutoff for determing iterative runs
      real    GL_i,GL_ii,GL_ir,GL_iir
      real    gapp(20), gapp45(20),gap_open

      integer np,L,niter,i,j
      integer invmap0(maxr),invmap(maxr)
      integer invmap_i (maxr),invmap_ii (maxr)
      integer invmap_ir(maxr),invmap_iir(maxr)
      integer best_init
      integer order_flag,score_flag,iflag,mode_flag,quick_flag
      integer mode,best_mode
      integer n_gapp,n_gapp45,i_gapp,id

      character*50 filename
      logical debug


      TMmax=0
      TM=0

      n_gapp=2
      gapp(1)=-0.6
      gapp(2)=0
      niter=30

      gapp45(1)=0
      n_gapp45=1

      TMcut1=-0.1
      TMcut2=0.2

      !coordinates of two structures
      do i=1,nseq1
         do j=1,3
            coor1(j,i)=xa(j,i,0)
         enddo
      enddo

      do i=1,nseq2
         do j=1,3
            coor2(j,i)=xa(j,i,1)
         enddo
      enddo

      best_init=-1
      best_mode=0
      iflag = 0
               

*111111111111111111111111111111111111111111111111111111111
*     initial alignment from gapless threading
**********************************************************
      call get_initial(iflag,nseq1,nseq2,ichainterm1(1),ichainterm2(1),
     &     id_chain1,id_chain2,invmap_i,invmap_ii,coor1,
     &     coor2,cont1,cont2,ncont1,ncont2,GL_i,GL_ii,nseq)

      call get_initial(iflag,nseq2,nseq1,ichainterm2(1),ichainterm1(1),
     &     id_chain2,id_chain1,invmap_ir,invmap_iir,coor2,
     &     coor1,cont2,cont1,ncont2,ncont1,GL_ir,GL_iir,nseq)

      !write(*,'(2(a,G))')'GL_i = ',GL_i,' GL_ii = ',GL_ii
      !write(*,'(2(a,G))')'GL_ir = ',GL_ir,' GL_iir = ',GL_iir

      if(GL_ii.lt.GL_iir)then
         call rev_invmap(invmap_iir,invmap_ii,nseq1,nseq2)
      endif

      if(GL_i.lt.GL_ir)then
         call rev_invmap(invmap_ir,invmap_i,nseq1,nseq2)
      endif

      ! save invmaps for 10
      do j=1,nseq2
         invmap_iir(j) = invmap_ii(j)
         invmap_ir(j)  = invmap_i(j)
      enddo

      TM_old=TMmax
      call iter_search(TMmax,TMcut1,score_flag,quick_flag,niter,
     &     gapp,n_gapp,invmap0,invmap_ii,nseq,nseq1,nseq2,12)
      if(TMmax-TM_old.gt.1e-6)then
         best_init=12
         TM_old=TMmax
      endif

      call iter_search(TMmax,TMcut2,score_flag,quick_flag,niter,
     &     gapp,n_gapp,invmap0,invmap_i,nseq,nseq1,nseq2,11)
      if(TMmax-TM_old.gt.1e-6)then
         best_init=11
         TM_old=TMmax
      endif

      !!! iterative non-sequential alignment
      if(order_flag.eq.0)then
         call iter_lsap(TMmax,score_flag,quick_flag,invmap_ii,invmap0,xa,nseq1,
     &        nseq2,nseq,ichainterm1,ichainterm2,3,mode_flag,best_mode,12)
         if(TMmax-TM_old.gt.1e-6)then
            best_init=12
            TM_old=TMmax
         endif

         call iter_lsap(TMmax,score_flag,quick_flag,invmap_i,invmap0,xa,nseq1,
     &        nseq2,nseq,ichainterm1,ichainterm2,3,mode_flag,best_mode,11)
         if(TMmax-TM_old.gt.1e-6)then
            best_init=11
            TM_old=TMmax
         endif
      endif

 100  continue 
*111111111111111111111111111100000000000000000000000000000000
*     iterative search first with TM-score,
*     then switch to IS-score
**********************************************************
      if(score_flag.eq.0)goto 222   !no need to repeat this for TM-score

****  search with TM-score
      BTMmax = 0

      call iter_search(BTMmax,TMcut1,iflag,quick_flag,3,gapp45,
     &     n_gapp45,invmap,invmap_iir,nseq,nseq1,nseq2,-1)

      call iter_search(BTMmax,TMcut2,iflag,quick_flag,3,gapp45,
     &     n_gapp45,invmap,invmap_ir,nseq,nseq1,nseq2,-1)

      !print *,'BTMmax =',BTMmax
      !call print_invmap(invmap,nseq2)

****  switch to IS-score
      call iter_search(TMmax,TMcut1,score_flag,quick_flag,3,
     &     gapp45,n_gapp45,invmap0,invmap,nseq,nseq1,nseq2,10)
      if(TMmax-TM_old.gt.1e-6)then
         best_init=10
         TM_old=TMmax
      endif
      
      if(order_flag.eq.0)then
         call iter_lsap(TMmax,score_flag,quick_flag,invmap,invmap0,xa,nseq1,
     &        nseq2,nseq,ichainterm1,ichainterm2,3,mode_flag,best_mode,10)
         if(TMmax-TM_old.gt.1e-6)then
            best_init=10
            TM_old=TMmax
         endif
      endif

 222  continue
*222222222222222222222222222222222222222222222222222222222222
*     get initial alignment from secondary structure alignment
**************************************************************
      call get_initial_from_ssc(isec,jsec,nseq1,nseq2,ichainterm1(1),
     &     ichainterm2(1),invmap0,invmap,xa,1.,0.01,nseq)

      call iter_search(TMmax,TMcut2,score_flag,quick_flag,niter,
     &     gapp45,n_gapp45,invmap0,invmap,nseq,nseq1,nseq2,2)
      if(TMmax-TM_old.gt.1e-6)then
         best_init=2
         TM_old=TMmax
      endif

      if(order_flag.eq.0)then
         call iter_lsap(TMmax,score_flag,quick_flag,invmap,invmap0,xa,nseq1,
     &        nseq2,nseq,ichainterm1,ichainterm2,3,mode_flag,best_mode,2)
         if(TMmax-TM_old.gt.1e-6)then
            best_init=2
            TM_old=TMmax
         endif
      endif


 333  continue
*3333333333333333333333333333333333333333333333333333333333333
*     get initial alignment of fragment superposition
**************************************************************
      call get_initial3(isec,jsec,nseq1,nseq2,ichainterm1(1),
     &     ichainterm2(1),id_chain1,id_chain2,invmap_i,coor1,
     &     coor2,cont1,cont2,ncont1,ncont2,GL_i,nseq,score_flag)

      !write(*,'(a,G)')'Init 3: GL_i = ',GL_i

      call iter_search(TMmax,TMcut2,score_flag,quick_flag,niter,
     &     gapp45,n_gapp45,invmap0,invmap_i,nseq,nseq1,nseq2,3)
      if(TMmax-TM_old.gt.1e-6)then
         best_init=3
         TM_old=TMmax
      endif

      !call print_invmap(invmap,nseq2)
      if(order_flag.eq.0)then
         call iter_lsap(TMmax,score_flag,quick_flag,invmap_i,invmap0,xa,nseq1,
     &        nseq2,nseq,ichainterm1,ichainterm2,5,mode_flag,best_mode,3)
         if(TMmax-TM_old.gt.1e-6)then
            best_init=3
            TM_old=TMmax
         endif
      endif

 444  continue
*444444444444444444444444444444444444444444444444444444444444
*     get initial alignment from invmap0+SS
*************************************************************
      call get_initial_from_ssc(isec,jsec,nseq1,nseq2,ichainterm1(1),
     &     ichainterm2(1),invmap0,invmap,xa,0.5,1.,nseq)


      call iter_search(TMmax,TMcut2,score_flag,quick_flag,niter,
     &     gapp45,n_gapp45,invmap0,invmap,nseq,nseq1,nseq2,4)
      if(TMmax-TM_old.gt.1e-6)then
         best_init=4
         TM_old=TMmax
      endif

      if(order_flag.eq.0)then
         call iter_lsap(TMmax,score_flag,quick_flag,invmap,invmap0,xa,nseq1,
     &        nseq2,nseq,ichainterm1,ichainterm2,3,mode_flag,best_mode,4)
         if(TMmax-TM_old.gt.1e-6)then
            best_init=4
            TM_old=TMmax
         endif
      endif


c^^^^^^^^^^^^^^^ best alignment invmap0(j) found ^^^^^^^^^^^^^^^^^^

 888  continue
!      write(*,'(a,a,I3,A8,F8.5,/)') 'Best alignment from', 
!     &     ' initial ',best_init, 'TMmax= ',TMmax



**888888888888888888888888888888888888888888888888888888888888888
*       initerative sequence-order-independent alignment
*****************************************************************

      if( order_flag.eq.1 ) goto 900  !skip sequence-order-independent iterative

      mode_flag=1
      do j=1,nseq2
         invmap(j)=invmap0(j)
      enddo
      debug=.false.

      call iter_lsap(TMmax,score_flag,quick_flag,invmap,invmap0,xa,nseq1,
     &     nseq2,nseq,ichainterm1,ichainterm2,20,mode_flag,best_mode,8)


 900  continue

      RETURN
      END
cccccc-------------- End of super_align ---------------------ccccc



cccccc======================================================ccccc
cccccc   Initerative searching for the best alignment
cccccc======================================================ccccc
      subroutine iter_search(TMmax,TMcut,score_flag,quick_flag,
     &  niter,gapp,n_gapp,invmap0,invmap_i,nseq,nseq1,nseq2,init)

      implicit none
      include  'pars.h'

      real     TM,TMmax,TM_old(10),rTMmax
      real     TMcut
      real     gapp(n_gapp),gap_open
      real     score(nseq,nseq),score0(nseq,nseq)
      real     diff,eps
      integer  niter,n_gapp,i_gapp
      integer  nseq1,nseq2,nseq
      integer  np,init
      integer  invmap(maxr), invmap0(maxr),invmap_i(maxr)
      integer  score_flag
      integer  i,j,id,ii,jj
      integer  mode,mode_flag,quick_flag

      logical  score_upd

      eps=1.E-6
      mode_flag=0

      do i=1,nseq2
         invmap(i)=invmap_i(i)
      enddo

      call get_score(nseq,invmap,TM,score0,score_flag,mode_flag,
     &     quick_flag,mode)                
      !write(*,'(a,i4,a,2F11.7)')'Initial ',init,' before DP',TM,TMmax

      if(TM - TMmax.gt.eps)then
         TMmax=TM
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
         do j=1,nseq2
            invmap_i(j)=invmap(j)
         enddo
      endif


      !for checking convergence and periodic cases, consider a period up to 5
      np=5
      do i=1,np
         TM_old(i)=-i  
      enddo

      !call print_invmap(invmap,nseq2)

      score_upd = .true.
      !Iterative dynamic programming
      if(TM-TMmax*TMcut.gt.eps.or.TMcut.lt.0)then
         DO 111 i_gapp=1,n_gapp	!different gap panalties
            GAP_OPEN=gapp(i_gapp) !gap open panalty, no gap extension penalty here
	    rTMmax=-1

            if(score_upd)then
               do ii=1,nseq1
                  do jj=1,nseq2
                     score(ii,jj)=score0(ii,jj)
                  enddo
               enddo
               score_upd=.false.
            endif

            do 222 id=1,niter      

               !dynamic programming. the alignment obtained is
               !saved in invmap(j), second seq -> first seq
               call DP(nseq,score,gap_open,nseq1,nseq2,invmap)


               !calculate TM-score/IS-score, score(i,j)
               call get_score(nseq,invmap,TM,score,score_flag,
     &              mode_flag,quick_flag,mode) 
               !write(*,'(F6.1,I4,2F11.7)')gapp(i_gapp),id,TM,TMmax


               if(TM-TMmax.gt.eps)then
                  TMmax=TM
                  do j=1,nseq2
                     invmap0(j)=invmap(j)
                  enddo
                  do j=1,nseq2
                     invmap_i(j)=invmap(j)
                  enddo
               endif

	       if(TM.gt.rTMmax)then
                  rTMmax=TM
                  do ii=1,nseq1
                     do jj=1,nseq2
                        score0(ii,jj)=score(ii,jj)
                     enddo
                  enddo
                  score_upd=.true.
     	       endif

               !check convergence and periodic cases
               do i=1,np
                  if(id.gt.i.and.abs(TM-TM_old(i)).lt.eps)goto 111
               enddo

               do i=np,2,-1
                  TM_old(i)=TM_old(i-1)
               enddo
               TM_old(1)=TM

 222        continue
 111     continue
      endif
      !write(*,'(a,i4,a,F11.7)')'Initial ',init,'  after DP: TMmax=',TMmax

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


cccccc======================================================ccccc
cccccc   Initerative non-sequential alignment search
cccccc======================================================ccccc
      subroutine iter_lsap(TMmax,score_flag,quick_flag,invmap,invmap0,xa,
     &     nseq1,nseq2,nseq,chterm1,chterm2,nround,mode_flag,best_mode,init)
      implicit none
      include 'pars.h'

      real score(nseq,nseq),new_score(nseq,nseq),z
      real TM,TMmax,TM_old(10)
      real xa(3,maxr,0:1)

      integer invmap(maxr),invmap0(maxr)
      integer invmap_i(maxr),invmap_ii(maxr),invmap1(maxr)
      integer chterm1(0:maxk),chterm2(0:maxk)
      integer nseq1,nseq2,nseq
      integer len11,len21,len12,len22,maxl
      integer score_flag,quick_flag,mode_flag
      integer nround,mode,best_mode
      integer id,np,i,j,init
      integer alnlen

      logical debug

      !call print_invmap(invmap,nseq2)

      debug=.false.
      mode_flag=1

      !Preparing scores for linear sum assignment
      call get_alnlen(invmap,nseq2,alnlen)
      if( alnlen.lt.3 ) return   ! must have at least three residues

      call get_score(nseq,invmap,TM,score,score_flag,mode_flag,
     &     quick_flag,mode)
      !write(*,'(i3,a,f11.7,a,f11.7)')init, ' Before iter_lsap: TM =',TM,',  TMamx = ',TMmax

      if(TM-TMmax.gt.1e-6)then
         TMmax=TM
         best_mode=mode
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
      endif

      !for checking convergence and periodic cases
      np=5 !period up to 5
      do i=1,np
         TM_old(i)=-i  
      enddo

      do 17 id=1,nround         !maximum interation

         !if(id.eq.1.and.init.eq.3) then
            !debug=.true.
            !call print_invmap(invmap,nseq2)
         !endif

ccc******** the first pair of chains
         len11=chterm1(1)
         len21=chterm2(1)
         call extract_score(nseq,score,nseq1,nseq2,1,1,len11,len21,new_score)

         maxl=max(len11,len21)
         call solvLSAP(nseq,new_score,invmap_i,len11,len21,maxl,debug)
         !debug=.false.


ccc******** the second pair of chains
         call extract_score(nseq,score,nseq1,nseq2,
     &        len11+1,len21+1,nseq1,nseq2,new_score)

         len12=nseq1-len11
         len22=nseq2-len21
         maxl=max(len12,len22)
         call solvLSAP(nseq,new_score,invmap_ii,len12,len22,maxl,debug)
         !debug=.false.

ccc******** now combines two optimal alignments
         call comb_invmap(invmap1,invmap_i,invmap_ii,len11,
     &        len21,len22,nseq2)

         call remove_distant_pairs(invmap1,invmap,nseq1,nseq2,xa)
         call get_alnlen(invmap,nseq2,alnlen)
         if( alnlen.lt.3 ) return  ! must have at least three residues

         call get_score(nseq,invmap,TM,score,score_flag,mode_flag,
     &        quick_flag,mode)
!         write(*,'(i3,a,i4,a,f11.7,a,f11.7,i4)')init,
!     &        ' round = ', id,', TM = ',TM,',  TMamx = ',TMmax, mode

         if(TM-TMmax.gt.1e-6)then
            TMmax=TM
            best_mode=mode
            do j=1,nseq2
               invmap0(j)=invmap(j)
            enddo
         endif

         !check convergence and periodic cases
         do i=1,np
            if(id.gt.i.and.abs(TM-TM_old(i)).lt.1e-6)goto 13
         enddo

         do i=np,2,-1
            TM_old(i)=TM_old(i-1)
         enddo
         TM_old(1)=TM

 17   continue
 13   continue

      return
      end


**************************************************************
*     get initial alignments from gapless threading
**************************************************************
      subroutine get_initial(score_flag,nseq1,nseq2,len1,len2,
     &     id_chain1,id_chain2,invmap_i,invmap_ii,coor1,coor2,
     &     cont1,cont2,ncont1,ncont2,GL_max,GL_maxA,nseq)
      implicit none
      include 'pars.h'

      real    coor1(3,maxr),coor2(3,maxr)
      real    xx,yy,zz,dd
      real    score(nseq,nseq),gap_open
      integer nseq1,nseq2,nseq   !length of full structures
      integer len1,len2     !length of first chains
      integer id_chain1(maxr),id_chain2(maxr)
      integer cont1(maxc,maxr),cont2(maxc,maxr)
      integer ncont1(maxr),ncont2(maxr)
      integer invmap(maxr),invmap_i(maxr),invmap_ii(maxr)

      common/d0/d0,anseq
      common/d0min/d0_min
      common/contcf/d_col,d_col2
      real d0,anseq,d0_min,d_col,d_col2
      real d01,d02
      real GL,GL_max,GL_maxA,GL_max_cut,aL

      double  precision r_1(3,maxr),r_2(3,maxr),r_3(3,maxr),w(maxr)
      double  precision u(3,3),t(3),rms,drms
      integer ier
      data w /maxr*1.0/
      
      real fcol                           !contact overlap factor
      integer imap1(maxr),imap2(maxr)     !aligned position to original position
      integer is_ali(maxr,maxr)  !original position to aligned position
      integer col
      integer score_flag
      integer alnlen, L,LL,ii,jj
      integer i,j,k,idel,n1,n2,ishift,n_jump

      ! initialize parameters and scores
      d01=d0+1.5
      if(d01.lt.d0_min)d01=d0_min
      d02=d01*d01
      n_jump=1


      GL_max=0
      GL_maxA=0
      GL_max_cut=0.95

      ! initialization 
      do i=1,nseq2
         invmap_i(i)=-1
         invmap_ii(i)=-1
      enddo

      aL=min(len1,len2)
      idel=aL/2.3             !minimum size of considered fragment
      if(idel.le.5)idel=5
      n1=-len2+idel
      n2=len1-idel

      call init_ali(imap1,imap2,nseq1,nseq2)

      do ishift=n1,n2,n_jump
         
         do k=1,nseq2
            invmap(k)=-1
         enddo

         L=0
         do 100 j=1,len2
            i=j+ishift
            if(i.ge.1.and.i.le.len1)then
               L=L+1
               invmap(j)=i
            endif
 100     continue

         call get_GL(GL,nseq1,nseq2,invmap,coor1,coor2,
     &        cont1,cont2,ncont1,ncont2,0)
         !print *, 'GL score =', GL

            
         if(GL.gt.GL_max*GL_max_cut)then  !this can lead to asymmetry
            LL=0
            do jj=1,len2
               ii=invmap(jj)
               if(ii.gt.0)then
                  LL=LL+1
                  r_1(1,LL)=coor1(1,ii)
                  r_1(2,LL)=coor1(2,ii)
                  r_1(3,LL)=coor1(3,ii)
                  r_2(1,LL)=coor2(1,jj)
                  r_2(2,LL)=coor2(2,jj)
                  r_2(3,LL)=coor2(3,jj)
               endif
            enddo
            call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2

             
ccc   is_ali(i,j) indicate whether i, j are in alignment dis<d_col
            do i=1,nseq1
               do j=1,nseq2
                  is_ali(i,j) = -1
               enddo
            enddo

            do i=1,nseq1
               xx=t(1)+u(1,1)*coor1(1,i)+u(1,2)*coor1(2,i)+u(1,3)*
     &              coor1(3,i)
               yy=t(2)+u(2,1)*coor1(1,i)+u(2,2)*coor1(2,i)+u(2,3)*
     &              coor1(3,i)
               zz=t(3)+u(3,1)*coor1(1,i)+u(3,2)*coor1(2,i)+u(3,3)*
     &              coor1(3,i)
               do j=1,nseq2
                  if(id_chain1(i).eq.id_chain2(j))then
                     dd=(xx-coor2(1,j))**2+(yy-coor2(2,j))**2
     &                    +(zz-coor2(3,j))**2
                     score(i,j)=1/(1+dd/d02)
                     if(dd.le.d_col2)then
                        is_ali(i,j)=1
                     endif
                  else
                     score(i,j)=ninf
                  endif
                  
               enddo
            enddo

ccc   calculate contact overlap and get IS score
            if(score_flag.eq.1) then
               do i=1,nseq1
                  do j=1,nseq2
                     if(id_chain1(i).eq.id_chain2(j))then
                        call get_fcol_4s(i,j,is_ali,fcol,col,
     &                       cont1,cont2,ncont1,ncont2)
                        score(i,j)=score(i,j)*(0.01+fcol)
                     endif
                  enddo
               enddo
            endif


*********extract alignement with score(i,j) *****************
            gap_open=0.
            call DP(nseq,score,gap_open,NSEQ1,NSEQ2,invmap)
            if(GL.gt.GL_max)then
               GL_max=GL
               do i=1,nseq2
                  invmap_i(i)=invmap(i)
               enddo
            endif

            call get_alnlen(invmap,nseq2,alnlen)
            !print *,'aln1=',alnlen,', idel=',idel
            if(alnlen.gt.idel)then

               call get_GL(GL,nseq1,nseq2,invmap,coor1,coor2,
     &              cont1,cont2,ncont1,ncont2,score_flag)
                  !print *,'GL=',GL,' GL_maxA=',GL_maxA
               if(GL.gt.GL_maxA)then
                  GL_maxA=GL
                  do i=1,nseq2
                     invmap_ii(i)=invmap(i)
                  enddo
               endif
            endif
         endif                  !GL_max condition ends
      enddo                     !threading loop ends

      return
      end


**************************************************************
*     get initial alignment from secondary structure alignment
**************************************************************
      subroutine get_initial_from_ssc(isec,jsec,nseq1,nseq2,
     &     len11,len21,invmap0,invmap,xa,w1,w2,nseq)
      implicit none
      include 'pars.h'

      integer invmap(maxr),invmap0(maxr)
      integer invmap_i(maxr),invmap_ii(maxr)
      integer isec(maxr),jsec(maxr)  !1->coil, 2->helix, 3->turn, 4->strand

      real    xa(3,maxr,0:1)
      real    score(nseq,nseq),score0(nseq,nseq),gap_open
      real    w1,w2   ! weight parameters for score assignment
      integer nseq1,nseq2,len11,len21,len12,len22,nseq
      integer i,j

      !get score(i,j) using RMSD martix
      call get_rms_score(nseq,score,invmap0,nseq1,nseq2,xa)


********** score matrix of chain 1s **************************
      do i=1,len11              !length of chain1, structure1
         do j=1,len21           !length of chain1, structure2
            if(isec(i).eq.jsec(j))then
               score0(i,j)=w1+w2*score(i,j) !avoid multiple optimal solutions
            else
               score0(i,j)=w2*score(i,j)
            endif
         enddo
      enddo

********** find initial alignment of chain1s ************
      gap_open=-1.0             !should be -1
      call DP(nseq,score0,gap_open,len11,len21,invmap_i)      !produce alignment for first chains


********** score matrix of chain 2s **************************
      len12 = nseq1 - len11  !length of chain2, structure1
      len22 = nseq2 - len21  !length of chain2, structure2

      do i=1,len12
         do j=1,len22
            if(isec(i+len11).eq.jsec(j+len21))then
               score0(i,j)=w1+w2*score(i+len11,j+len21) !w2 avoid multiple optimal solutions
            else
               score0(i,j)=w2*score(i+len11,j+len21)
            endif
         enddo
      enddo

********** find initial alignment of chain 2s ************
      call DP(nseq,score0,gap_open,len12,len22,invmap_ii)   !produce alignment for second chains

      !call print_invmap(invmap_i,len21)
      !call print_invmap(invmap_ii,len22)

********** piece together two alignments ************
      do j=1,nseq2
         invmap(j)=-1
      enddo

      do j=1,len21
         i=invmap_i(j)
         if(i.gt.0) invmap(j)=i
      enddo
      do j=1,len22
         i=invmap_ii(j)
         if(i.gt.0) invmap(j+len21)=i+len11
      enddo

      !call print_invmap(invmap,nseq2)

      return
      end




**************************************************************
*    Fragment-based superposition.                           *
**************************************************************
      subroutine get_initial3(sec1,sec2,nseq1,nseq2,len1,len2,
     &     id_chain1,id_chain2,invmap_i,coor1,coor2,
     &     cont1,cont2,ncont1,ncont2,GLmaxA,nseq,score_flag)
      implicit none
      include 'pars.h'

      integer id_chain1(maxr),id_chain2(maxr)
      integer cont1(maxc,maxr),cont2(maxc,maxr)     !contact list
      integer ncont1(maxr),ncont2(maxr)      !number of contacts
      integer sec1(maxr),sec2(maxr)  !1->coil, 2->helix, 3->turn, 4->strand
      integer invmap(maxr),invmap_i(maxr)
      integer nseq1,nseq2,nseq
      integer len1,len2,minlen,maxlen
      integer score_flag

      real    coor1(3,maxr),coor2(3,maxr)
      real    score(nseq,nseq),gap_open

      common/d0/d0,anseq
      common/d0min/d0_min
      common/contcf/d_col,d_col2
      real    d0,anseq,d0_min,d01,d02,dd
      real    xx,yy,zz
      real    d_col,d_col2
      real    fcol               !contact overlap factor
      integer is_ali(maxr,maxr)  !original position to aligned position
      integer col

      double precision r_1(3,maxr),r_2(3,maxr),r_3(3,maxr),w(maxr)
      double precision u(3,3),t(3),rms
      integer ier
      data w /maxr*1.0/

      integer n_frag,m1,m2,alnlen
      integer n_jump1,n_jump2
      integer i,j,ii,jj,k
      integer score_flag_i
      integer secmatch,n_frag2

      real GL,GL_cf,GLmaxA,GLmaxA_cf

     
***** setting parameters ************************************

      maxlen=max(len1,len2)
      minlen=min(len1,len2)

      n_frag=minlen/20   !length of the sliding window
      if(n_frag.lt.5)n_frag=5

      n_jump1=len1/20 - 2 
      n_jump2=len2/20 - 2
      if(n_jump1.lt.1)n_jump1=1
      if(n_jump2.lt.1)n_jump2=1
      if(n_jump1.gt.n_frag)n_jump1=n_frag
      if(n_jump2.gt.n_frag)n_jump2=n_frag

      m1=len1-n_frag+1
      m2=len2-n_frag+1

      n_frag2=n_frag/2
      if(n_frag2.lt.2)n_frag2=2


      d01=d0+1.5
      if(d01.lt.d0_min)d01=d0_min
      d02=d01*d01

      GLmaxA=-1
      GLmaxA_cf=min(nseq1,nseq2)*0.9
      GL_cf=n_frag/2.

      score_flag_i=0  !TM-score

      !initialization 
      do i=1,nseq2
         invmap_i(i)=-1
      enddo

      !main loop, fragment vs fragment
      do 10 ii=1,m1,n_jump1
         do 20 jj=1,m2,n_jump2
            secmatch=0
            do j=1,nseq2
               invmap(j)=-1
            enddo

            do k=1,n_frag
               i=ii+k-1
               j=jj+k-1
               invmap(j)=i
               if(sec1(i).eq.sec2(j))then
                  secmatch=secmatch+1
               endif
            enddo

            if(secmatch.le.n_frag2) goto 20

            !!! first identify similar fragements
            call get_GL(GL,nseq1,nseq2,invmap,coor1,coor2,
     &           cont1,cont2,ncont1,ncont2,score_flag_i)

            if(GL.le.GL_cf)goto 20

            !!! use good fragment alignment to rotate whole structures
            k=0
            do 200 j=1,nseq2
               i=invmap(j)
               if(i.le.0)goto 200
               k=k+1
               r_1(1,k)=coor1(1,i)
               r_1(2,k)=coor1(2,i)
               r_1(3,k)=coor1(3,i)
               r_2(1,k)=coor2(1,j)
               r_2(2,k)=coor2(2,j)
               r_2(3,k)=coor2(3,j)
 200        continue

ccc   is_ali(i,j) indicate whether i, j are in alignment dis<d_col
            do i=1,nseq1
               do j=1,nseq2
                  is_ali(i,j) = -1
               enddo
            enddo

*********superpose the two structures and rotate it *****************
            call u3b(w,r_1,r_2,k,1,rms,u,t,ier) !u rotate r_1 to r_2
            do i=1,nseq1
              xx=t(1)+u(1,1)*coor1(1,i)+u(1,2)*coor1(2,i)+u(1,3)*coor1(3,i)
              yy=t(2)+u(2,1)*coor1(1,i)+u(2,2)*coor1(2,i)+u(2,3)*coor1(3,i)
              zz=t(3)+u(3,1)*coor1(1,i)+u(3,2)*coor1(2,i)+u(3,3)*coor1(3,i)
              do j=1,nseq2
                 if(id_chain1(i).eq.id_chain2(j))then
                    dd=(xx-coor2(1,j))**2+(yy-coor2(2,j))**2
     &                   +(zz-coor2(3,j))**2
                    score(i,j)=1/(1+dd/d02)
                    if(dd.le.d_col2)then
                       is_ali(i,j)=1
                    endif
                 else
                    score(i,j)=ninf
                 endif
              enddo
           enddo
                  
ccc   calculate contact overlap and get IS score
           if(score_flag.eq.1) then
              do i=1,nseq1
                 do j=1,nseq2
                    if(id_chain1(i).eq.id_chain2(j))then
                       call get_fcol_4s(i,j,is_ali,fcol,col,
     &                      cont1,cont2,ncont1,ncont2)
                       score(i,j)=score(i,j)*(0.01+fcol)
                    endif
                 enddo
              enddo
           endif
*********extract alignement with score(i,j) *****************

           gap_open=0
           call DP(nseq,score,gap_open,NSEQ1,NSEQ2,invmap)

           call get_GL(GL,nseq1,nseq2,invmap,coor1,coor2,
     &          cont1,cont2,ncont1,ncont2,score_flag)

           if(GL.gt.GLmaxA)then
              GLmaxA=GL
              do j=1,nseq2
                 invmap_i(j)=invmap(j)
              enddo
              if(GLmaxA.gt.GLmaxA_cf)return  !reduce scans for highly similar structures
           endif

 20     continue
 10   continue
      
      return
      end







****************************************************************
*     calculate TM-score and score matrix for DP and LSAP 
****************************************************************
      subroutine get_score(nseq,invmap,TM,score,score_flag,
     &     mode_flag,quick_flag,mode)
      implicit none
      include 'pars.h'

      common/length/nseq1,nseq2
      COMMON/BACKBONE/XA(3,maxr,0:1)
      common/d0/d0,anseq
      common/contcf/d_col,d_col2

      common/id2chain/id_chain1(maxr),id_chain2(maxr)
      common/contlst/cont1,cont2
      common/contnum/ncont1(maxr),ncont2(maxr)  !number of contacts
      integer cont1(maxc,maxr),cont2(maxc,maxr) !contact list
      integer id_chain1,id_chain2,ncont1,ncont2

      integer nseq1,nseq2,nseq,n_al
      integer invmap(maxr),invmap_r(maxr)

      real    xa,d0,anseq,d_col,d_col2
      real    xtm1(maxr),ytm1(maxr),ztm1(maxr)
      real    xtm2(maxr),ytm2(maxr),ztm2(maxr)
      real    xtm1r(maxr),ytm1r(maxr),ztm1r(maxr)
      real    xtm2r(maxr),ytm2r(maxr),ztm2r(maxr)
      real    xx,yy,zz,dd
      real    score(nseq,nseq)
      real    TM, TMr, Rcomm, diff,diff1,diff2

      real    fcol              !contact overlap factor
      real    d0_input,d0_search
      integer is_ali(maxr,maxr) !original position to aligned position
      integer imap1(maxr),imap2(maxr)
      integer col,ncol,ncolr,Lcomm
      integer score_flag,mode_flag,mode
      integer quick_flag,isearch

ccc   RMSD:
      double precision r_1(3,maxr),r_2(3,maxr),r_1r(3,maxr),w(maxr)
      double precision u(3,3),t(3),rms,drms !armsd is real
      integer ier
      data w /maxr*1.0/
ccc   
      integer i,j

      mode=0
      isearch=2  !extensive search
      if(quick_flag.ge.1) isearch=1


c     calculate RMSD between aligned structures and rotate the structures -->
      n_al=0
      do j=1,NSEQ2
         i=invmap(j)            !j aligned to i
         if(i.gt.0)then
            n_al=n_al+1
ccc   for TM-score:
            xtm1(n_al)=xa(1,i,0) !for TM-score
            ytm1(n_al)=xa(2,i,0)
            ztm1(n_al)=xa(3,i,0)
            xtm2(n_al)=xa(1,j,1)
            ytm2(n_al)=xa(2,j,1)
            ztm2(n_al)=xa(3,j,1)
ccc   for rotation matrix:
            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
            imap1(n_al)=i       !aligned position to original index
            imap2(n_al)=j
         endif
      enddo

      TM=-1
      if(n_al.lt.3)return  !must have at least three residues


***   calculate TM-score for the given alignment----------->
      d0_input=d0
      d0_search=d0
      if (d0_search.gt.8.0)d0_search=8.0
      if (d0_search.lt.1.0)d0_search=1.0

      call TMsearch(d0_input,d0_search,n_al,xtm1,ytm1,ztm1,imap1,
     &     n_al,xtm2,ytm2,ztm2,imap2,TM,Rcomm,Lcomm,isearch,
     &     score_flag,ncol)     !simplified search engine
      TM=TM*n_al/anseq          !TM-score


ccccc invmap_r maps sequence 1 to sequence 2, and use it to
ccccc calculate TM-score again. note that TM-score calculation is
ccccc residue order dependent. The condition block below ensures that
ccccc TM-score is symmetrical during the non-squential alignment search
      if(mode_flag.eq.1)then
         do i=1,nseq1
            invmap_r(i)=-1   !seq1 to seq2, reverse of invmap
         enddo

         do j=1,nseq2
            i=invmap(j)
            if(i.gt.0) invmap_r(i)=j
         enddo

         n_al=0
         do i=1,nseq1
            j=invmap_r(i)         !j aligned to i
            if(j.gt.0)then
               n_al=n_al+1
               xtm1r(n_al)=xa(1,i,0) !for TM-score
               ytm1r(n_al)=xa(2,i,0)
               ztm1r(n_al)=xa(3,i,0)
               xtm2r(n_al)=xa(1,j,1)
               ytm2r(n_al)=xa(2,j,1)
               ztm2r(n_al)=xa(3,j,1)

               r_1r(1,n_al)=xa(1,i,0)
               r_1r(2,n_al)=xa(2,i,0)
               r_1r(3,n_al)=xa(3,i,0)

               imap1(n_al)=i    !aligned position to original index
               imap2(n_al)=j
            endif
         enddo

         call TMsearch(d0_input,d0_search,n_al,xtm1r,ytm1r,ztm1r,
     &        imap1,n_al,xtm2r,ytm2r,ztm2r,imap2,TMr,Rcomm,Lcomm,
     &        isearch,score_flag,ncolr)  !simplified search engine
         TMr=TMr*n_al/anseq       !TM-score

         !!! ensure the alignment is same as that from two structures in revsersed order
         diff=TMr-TM
         diff1=abs(xa(1,1,0)-xa(1,2,0)) ! arbitrarily break symmetry
         diff2=abs(xa(1,1,1)-xa(1,2,1))
         if(diff.gt.1e-6.or.(abs(diff).le.1e-6.and.diff1.gt.diff2))then
            TM=TMr
            ncol=ncolr
            mode=1
            do i=1,n_al
               xtm1(i)=xtm1r(i)
               ytm1(i)=ytm1r(i)
               ztm1(i)=ztm1r(i)

               r_1(1,i)=r_1r(1,i)
               r_1(2,i)=r_1r(2,i)
               r_1(3,i)=r_1r(3,i)
            enddo
         endif
      endif




ccc   is_ali(i,j) indicate whether i, j are in alignment dis<d_col
      do i=1,nseq1
         do j=1,nseq2
            is_ali(i,j) = -1
         enddo
      enddo


***   calculate score matrix score(i,j)------------------>
      do i=1,n_al
         r_2(1,i)=xtm1(i)
         r_2(2,i)=ytm1(i)
         r_2(3,i)=ztm1(i)
      enddo

      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) !u rotate r_1 to r_2

      do i=1,nseq1
         xx=t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+u(1,3)*xa(3,i,0)
         yy=t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+u(2,3)*xa(3,i,0)
         zz=t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+u(3,3)*xa(3,i,0)
         do j=1,nseq2
            if(id_chain1(i).eq.id_chain2(j))then
               dd=(xx-xa(1,j,1))**2+(yy-xa(2,j,1))**2+(zz-xa(3,j,1))**2
               score(i,j)=1/(1+dd/d0**2)

               if(dd.le.d_col2)then
                  is_ali(i,j)=1
               endif
            else
               score(i,j)=ninf
            endif
         enddo
      enddo

ccc   calculate contact overlap and get IS score
      ncol=0  ! total number of contact overlaps
      do i=1,nseq1
         do j=1,nseq2
            if(score_flag.eq.1.and.id_chain1(i).eq.id_chain2(j))then
               call get_fcol_4s(i,j,is_ali,fcol,col,
     &              cont1,cont2,ncont1,ncont2)
               ncol=ncol+col
               score(i,j)=score(i,j)*(fcol+0.01)
            endif
            !print *,i,j,score(i,j),fcol
         enddo
      enddo
      ncol=ncol/2  !each contact counted twice
      
      !print *,'get_score is done'
c^^^^^^^^^^^^^^^^ score(i,j) done ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end


****************************************************************
*     remove Ca pairs with distance larger than d8 from invmap 
****************************************************************
      subroutine remove_distant_pairs(invmap,invmap_old,nseq1,nseq2,xa)
      implicit none
      include 'pars.h'

      integer nseq1,nseq2,nseq
      integer invmap(maxr),invmap_old(maxr)

      common/d8/d8
      real    d8

      real    xtm1(maxr),ytm1(maxr),ztm1(maxr)
      real    xtm2(maxr),ytm2(maxr),ztm2(maxr)
      real    xx,yy,zz
      real    XA(3,maxr,0:1)
      real    dd

ccc   RMSD:
      double precision r_1(3,maxr),r_2(3,maxr),r_3(3,maxr),w(maxr)
      double precision u(3,3),t(3),rms,drms !armsd is real
      integer ier
      data w /maxr*1.0/
ccc   
      integer i,j,n_al

c     calculate RMSD between aligned structures and rotate the structures -->
      n_al=0
      do j=1,nseq2
         i=invmap_old(j)            !j aligned to i
         if(i.gt.0)then
            n_al=n_al+1
ccc   for rotation matrix:
            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
            r_2(1,n_al)=xa(1,j,1)
            r_2(2,n_al)=xa(2,j,1)
            r_2(3,n_al)=xa(3,j,1)
         endif
      enddo
***   calculate score matrix score(i,j)------------------>
      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) !u rotate r_1 to r_2

      do j=1,nseq2
         invmap_old(j) = -1
      enddo

      do 10 j=1,nseq2
         i=invmap(j)
         if(i.le.0) goto 10
         xx=t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+u(1,3)*xa(3,i,0)
         yy=t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+u(2,3)*xa(3,i,0)
         zz=t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+u(3,3)*xa(3,i,0)

         dd=sqrt((xx-xa(1,j,1))**2+(yy-xa(2,j,1))**2+(zz-xa(3,j,1))**2)
         if(dd.le.d8) then
            invmap_old(j) = i
         endif
 10   continue

c^^^^^^^^^^^^^^^^ score(i,j) done ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end




****************************************************************
*     with invmap(i) calculate score(i,j) using RMSD rotation
****************************************************************
      subroutine get_rms_score(nseq,score,invmap,nseq1,nseq2,xa)
      implicit none
      include 'pars.h'

      integer nseq1,nseq2,nseq
      integer invmap(maxr)
      integer id_chain1(maxr),id_chain2(maxr)

      common/d0/d0,anseq
      common/d0min/d0_min
      real    d0,anseq,d0_min,d01,d02

      real    xtm1(maxr),ytm1(maxr),ztm1(maxr)
      real    xtm2(maxr),ytm2(maxr),ztm2(maxr)
      real    xx,yy,zz
      real    XA(3,maxr,0:1)
      real    score(nseq,nseq)
      real    dd

ccc   RMSD:
      double precision r_1(3,maxr),r_2(3,maxr),r_3(3,maxr),w(maxr)
      double precision u(3,3),t(3),rms,drms !armsd is real
      integer ier
      data w /maxr*1.0/
ccc
      integer i,j,n_al

c     calculate RMSD between aligned structures and rotate the structures -->
      n_al=0
      do j=1,NSEQ2
         i=invmap(j)            !j aligned to i
         if(i.gt.0)then
            n_al=n_al+1
ccc   for rotation matrix:
            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
            r_2(1,n_al)=xa(1,j,1)
            r_2(2,n_al)=xa(2,j,1)
            r_2(3,n_al)=xa(3,j,1)
         endif
      enddo
***   calculate score matrix score(i,j)------------------>
      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) !u rotate r_1 to r_2
      d01=d0+1.5
      if(d01.lt.d0_min)d01=d0_min
      d02=d01*d01
      do i=1,nseq1
         xx=t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+u(1,3)*xa(3,i,0)
         yy=t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+u(2,3)*xa(3,i,0)
         zz=t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+u(3,3)*xa(3,i,0)
         do j=1,nseq2
            if(id_chain1(i).eq.id_chain2(j))then
               dd=(xx-xa(1,j,1))**2+(yy-xa(2,j,1))**2+(zz-xa(3,j,1))**2
               score(i,j)=1/(1+dd/d02)
            else
               score(i,j)=ninf
            endif
         enddo
      enddo

c^^^^^^^^^^^^^^^^ score(i,j) done ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end





********************************************************************
*     dynamic programming for alignment.
*     input: score(i,j), and gap_open
*     output: invmap(j)
********************************************************************
      subroutine DP(n,score,gap_open,len_seq1,len_seq2,invmap)
      implicit none
      include 'pars.h'

      logical*1 dir
      real    score(n,n),gap_open
      integer n,invmap(maxr),len_seq1,len_seq2
      integer i,j

      dimension dir(0:n,0:n),val(0:n,0:n)
      real*8   h,v,d,val
      


***   initialize the matrix with no end gap penalty:
      val(0,0)=0
      do i=1,len_seq1
        dir(i,0)=.false.
        val(i,0)=0
      enddo
      do j=1,len_seq2
        dir(0,j)=.false.
        val(0,j)=0
        invmap(j)=-1
      enddo

***   calculate the F matrix and path:
      do j=1,len_seq2
        do i=1,len_seq1
          d=val(i-1,j-1)+score(i,j)
          h=val(i-1,j)
          if(dir(i-1,j))h=h+gap_open
          v=val(i,j-1)
          if(dir(i,j-1))v=v+gap_open
          
          if((d.ge.h).and.(d.ge.v)) then
            dir(i,j)=.true.
            val(i,j)=d
          else
            dir(i,j)=.false.
            if(v.gt.h)then
              val(i,j)=v
            else
              val(i,j)=h
            end if
          endif
        enddo
      enddo
      
      !print *,'optimal dp cost is',val(len_seq1,len_seq2)
***   extract the alignment:
      i=len_seq1
      j=len_seq2
      do while((i.gt.0).and.(j.gt.0))
        if(dir(i,j))then
          invmap(j)=i
          i=i-1
          j=j-1
        else
          h=val(i-1,j)
          if(dir(i-1,j))h=h+gap_open
          v=val(i,j-1)
          if(dir(i,j-1))v=v+gap_open
          if(v.gt.h) then
            j=j-1
          else
            i=i-1
          endif
        endif
      enddo
      

      !print *,'optimal dp cost: ',val(len_seq1,len_seq2)
c^^^^^^^^^^^^^^^dynamical programming done ^^^^^^^^^^^^^^^^^^^
      return
      end


**************************************************************
*     calculate secondary structure from Calpha atoms
**************************************************************
      subroutine calcSEC(strid,sec,nseq)
      include 'pars.h'
      integer nseq,sec(maxr)
      integer strid

********** assign secondary structures ***************
c     1->coil, 2->helix, 3->turn, 4->strand
      do i=1,nseq
         sec(i)=1
         j1=i-2
         j2=i-1
         j3=i
         j4=i+1
         j5=i+2
         if(j1.ge.1.and.j5.le.nseq)then
            dis13=diszy(strid,j1,j3)
            dis14=diszy(strid,j1,j4)
            dis15=diszy(strid,j1,j5)
            dis24=diszy(strid,j2,j4)
            dis25=diszy(strid,j2,j5)
            dis35=diszy(strid,j3,j5)
            sec(i)=make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
         endif
      enddo

      call smooth(sec,nseq)               !smooth the assignment

      return
      end


**************************************************************
*     smooth the secondary structure assignment
**************************************************************
      subroutine smooth(isec,nseq1)
      include 'pars.h'

      integer isec(maxr),nseq1

***   smooth single -------------->
***   --x-- => -----
      do i=1,nseq1
         if(isec(i).eq.2.or.isec(i).eq.4)then
            j=isec(i)
            if(isec(i-2).ne.j)then
               if(isec(i-1).ne.j)then
                  if(isec(i+1).ne.j)then
                     if(isec(i+2).ne.j)then  ! Mu: i+2 instead of i+1
                        isec(i)=1
                     endif
                  endif
               endif
            endif
         endif
      enddo


***   smooth double -------------->
***   --xx-- => ------
      do i=1,nseq1
         if(isec(i).ne.2)then
         if(isec(i+1).ne.2)then
         if(isec(i+2).eq.2)then
         if(isec(i+3).eq.2)then
         if(isec(i+4).ne.2)then
         if(isec(i+5).ne.2)then
            isec(i+2)=1
            isec(i+3)=1
         endif
         endif
         endif
         endif
         endif
         endif

         if(isec(i).ne.4)then
         if(isec(i+1).ne.4)then
         if(isec(i+2).eq.4)then
         if(isec(i+3).eq.4)then
         if(isec(i+4).ne.4)then
         if(isec(i+5).ne.4)then
            isec(i+2)=1
            isec(i+3)=1
         endif
         endif
         endif
         endif
         endif
         endif
      enddo


***   connect -------------->
***   x-x => xxx
      do i=1,nseq1
         if(isec(i).eq.2)then
         if(isec(i+1).ne.2)then
         if(isec(i+2).eq.2)then
            isec(i+1)=2
         endif
         endif
         endif

         if(isec(i).eq.4)then
         if(isec(i+1).ne.4)then
         if(isec(i+2).eq.4)then
            isec(i+1)=4
         endif
         endif
         endif
      enddo

      return
      end

*************************************************************
*     assign secondary structure:
*************************************************************
      function diszy(i,i1,i2)
      include 'pars.h'

      COMMON/BACKBONE/XA(3,maxr,0:1)
      diszy=sqrt((xa(1,i1,i)-xa(1,i2,i))**2
     &     +(xa(2,i1,i)-xa(2,i2,i))**2
     &     +(xa(3,i1,i)-xa(3,i2,i))**2)
      return
      end

*************************************************************
*     assign secondary structure:
*************************************************************
      function make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
      make_sec=1
      delta=2.1
      if(abs(dis15-6.37).lt.delta)then
         if(abs(dis14-5.18).lt.delta)then
            if(abs(dis25-5.18).lt.delta)then
               if(abs(dis13-5.45).lt.delta)then
                  if(abs(dis24-5.45).lt.delta)then
                     if(abs(dis35-5.45).lt.delta)then
                        make_sec=2 !helix
                        return
                     endif
                  endif
               endif
            endif
         endif
      endif
      delta=1.42
      if(abs(dis15-13).lt.delta)then
         if(abs(dis14-10.4).lt.delta)then
            if(abs(dis25-10.4).lt.delta)then
               if(abs(dis13-6.1).lt.delta)then
                  if(abs(dis24-6.1).lt.delta)then
                     if(abs(dis35-6.1).lt.delta)then
                        make_sec=4 !strand
                        return
                     endif
                  endif
               endif
            endif
         endif
      endif
      if(dis15.lt.8)then
         make_sec=3
      endif

      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Read a contact file                                              c
c     The format of file is                                            c
c     Res  ResidueIndex  NumContacts  Contact Residues                 c
c     e.g.
c     RES    74   3   49   50   51                                     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readContFile(contfile,ncont,cont)
      include 'pars.h'

      character*500 contfile
      character*1000 line
      character*5,head
      integer ncont(maxr),cont(maxc,maxr)


****** Initializing contact arrays ********
      do i=1,maxr
         ncont(i)=0
         do j=1,maxc
            cont(j,i)=0
         enddo
      enddo

      open(unit=15,file=contfile,status='old')
      do while (.true.)
         read(15,1001,end=100) line
         if(line(1:5).eq.'RES  ') then
            read(line,*) head, ires, ncont(ires),
     &           (cont(i,ires),i=1,ncont(ires))
         endif
      enddo
 100  continue
 1001 format(A1000)

      close(15)
      return
      end


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Calculate residue-residue contact at interface                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine calcContact(xa,strid,chterm,cont,ncont,cf)
      implicit none
      include 'pars.h'

      real xa(3,maxr,0:1)
      real cf,cf2
      real dis2

      integer cont(maxc,maxr),ncont(maxr),chterm(0:maxk)
      integer i,j,strid

      cf2=cf*cf

****** Initializing contact arrays ********
      do i=1,maxr
         ncont(i)=0
         do j=1,maxc
            cont(j,i)=0
         enddo
      enddo

****** calculate contacts ******
      do i=1,chterm(1)
         do j=chterm(1)+1,chterm(2)
            dis2=(xa(1,i,strid)-xa(1,j,strid))**2+
     &           (xa(2,i,strid)-xa(2,j,strid))**2+
     &           (xa(3,i,strid)-xa(3,j,strid))**2
            if(dis2.lt.cf2)then
               ncont(i)=ncont(i)+1
               ncont(j)=ncont(j)+1
               cont(ncont(i),i)=j
               cont(ncont(j),j)=i
            endif
         enddo
      enddo

      return
      end



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Read a PDB file.   Please check the PDB coordinate format at:    c
c     http://www.wwpdb.org/documentation/format32/sect9.html           c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readPDBFile(pdbfile,coor,id,nchain,ichterm,
     &     chname,seq,respdbnm,respdbid,sec,length,cont_flag,cont,ncont)
      implicit none
      include 'pars.h'

      character*100 pdbfile, line
      character*5  respdbid(maxr)
      character*3  aa3(maxr),respdbnm(maxr)
      character*1  chname(maxk),seq(maxr)
      character*1 altloc,chain,chtmp
      character*3 resnam,chtmp3
      character*5 atmnam,resid
      character*6 head

      real   coor(3,maxr,0:1)
      real   occupancy(maxr)  !occupancy
      real   ca_ca_cf
      real   x,y,z,occ,temp

      integer  nchain, ichterm(0:maxk)
      integer  length     !total number of residues
      integer  cont_flag,id
      integer  sec(maxr),cont(maxc,maxr),ncont(maxr)
      integer  iatom,ires
      integer  i,j,k,nc,counter,imap(maxr)

      character*3 aa(-1:19)
      data aa/ 'UNK','GLY','ALA','SER','CYS',
     &     'VAL','THR','ILE','PRO','MET',
     &     'ASP','ASN','LEU','LYS','GLU',
     &     'GLN','ARG','HIS','PHE','TYR',
     &     'TRP'/

      character*1 slc(-1:19)
      data slc/'X','G','A','S','C',
     &     'V','T','I','P','M',
     &     'D','N','L','K','E',
     &     'Q','R','H','F','Y',
     &     'W'/

      open(unit=8,file=pdbfile,status='old')

      i=0
      nchain=0
      chain=' '
      ichterm(0)=0
      do while (.true.)
         read(8,9001,end=100) line

         if(line(1:3).eq.'TER')then
            nchain=nchain+1     ! count number of chains
            ichterm(nchain)=i   ! index of the last residue of a chain
            if(chain.eq.' ') chain='_'
            chname(nchain)=chain
         endif

         if(line(1:4).eq.'ATOM')then
            read(line,9000) head,iatom,chtmp,atmnam,altloc,resnam,
     &           chtmp,chain,resid,chtmp3,x,y,z,occ,temp

            if(atmnam.ne.'CA  '.and.atmnam.ne.' CA '
     &           .and.atmnam.ne.'  CA') goto 90

            if(altloc.ne.' '.and.altloc.ne.'A'.and.altloc.ne.'1')
     &           goto 90

            i=i+1
            coor(1,i,id)=x
            coor(2,i,id)=y
            coor(3,i,id)=z
            respdbnm(i)=resnam
            respdbid(i)=resid
            occupancy(i)=occ
            sec(i)=nint(occ)
            
            do j=-1,19
               if(resnam.eq.aa(j))then
                  seq(i)=slc(j)
                  goto 21
               endif
            enddo
            seq(i)=slc(-1)
 21         continue

            !print "(I5,A2,2A6)",i,seq(i),respdbnm(i),respdbid(i)
            if(i.ge.maxr)then
               print *, "Warning: total number of residues beyond ",maxr
               goto 100
            endif
         endif

 90      continue
      enddo
 100  continue
 9000 format(A6,I5,A1,A4,A1,A3,2A1,A5,A3,3F8.3,2F6.2)
 9001 format(A100)

      close(8)
      length=i

      if(cont_flag.eq.1)then
         call calcSEC(id,sec,length)  ! secondary structure

         ca_ca_cf=8.5    ! distace cutoff 8.5 A
         call calcContact(coor,id,ichterm,cont,ncont,ca_ca_cf)  !interface contacts

         ! select only interface residues with at least one contact
         nc=0
         do i=1,length
            imap(i)=0
            if(ncont(i).gt.0)then
               nc=nc+1
               coor(1,nc,id)=coor(1,i,id)
               coor(2,nc,id)=coor(2,i,id)
               coor(3,nc,id)=coor(3,i,id)
               respdbnm(nc)=respdbnm(i)
               respdbid(nc)=respdbid(i)
               sec(nc)=sec(i)
               seq(nc)=seq(i)
               imap(i)=nc
            endif
            do j=1,nchain
               if(ichterm(j).eq.i)ichterm(j)=nc
            enddo
         enddo

         ! re-index contacts
         do 10 i=1,length
            if(imap(i).eq.0) goto 10
            counter=0
            do k=1,ncont(i)
               j=cont(k,i)
               if(imap(j).gt.0)then
                  counter=counter+1
                  cont(counter,imap(i))=imap(j)
               endif
            enddo
            ncont(imap(i))=counter
 10      continue

         length=nc
      endif


      return
      end





cccc=======================================================================cccc
cccc  Swap two chains in a PDB file                                        cccc
cccc  This may be extended to generate arbitrary combination               cccc
cccc  for mulitiple chains                                                 cccc
cccc=======================================================================cccc
      subroutine swap_chains(nchain,chterm,chname,seq,coor,pdbid,
     &     respdbnm,respdbid,sec,id_chain,ncont,cont)

      include 'pars.h'

      character*5 respdbid(maxr), respdbid_tmp(maxr)
      character*3 respdbnm(maxr), respdbnm_tmp(maxr)
      character   chname(maxk),chname_tmp(maxk)
      character   seq(maxr), seq_tmp(maxr)

      real coor(3,maxr,0:1),coor_tmp(3,maxr,0:1)

      integer nchain, chterm(0:maxk), chterm_tmp(0:maxk)
      integer id_chain(maxr),sec(maxr), sec_tmp(maxr)
      integer ncont(maxr),ncont_tmp(maxr)
      integer cont(maxc,maxr), cont_tmp(maxc,maxr)
      integer pdbid, pos

      if(nchain.ne.2) return


      !store all variables to temp places
      do i=1,chterm(nchain)
         respdbid_tmp(i) = respdbid(i)
         respdbnm_tmp(i) = respdbnm(i)
         seq_tmp(i) = seq(i)

         !print "(A4,I4,1X,A1,A6)",'i= ',i,seq(i),respdbnm(i)
         do j=1,3
            coor_tmp(j,i,pdbid) = coor(j,i,pdbid)
         enddo

         sec_tmp(i)=sec(i)
         ncont_tmp(i)=ncont(i)
         do j=1,ncont(i)
            cont_tmp(j,i)=cont(j,i)
         enddo
      enddo

      do i=1,nchain
         chname_tmp(i)=chname(i)
         chterm_tmp(i)=chterm(i)
      enddo

      !now swapping
      pos=0
      do i=2,1,-1
         do j=chterm_tmp(i-1)+1,chterm_tmp(i)
            pos=pos+1
            respdbid(pos)=respdbid_tmp(j)
            respdbnm(pos)=respdbnm_tmp(j)
            seq(pos)=seq_tmp(j)

            do k=1,3
               coor(k,pos,pdbid)=coor_tmp(k,j,pdbid)
            enddo

            sec(pos)=sec_tmp(j)
            ncont(pos)=ncont_tmp(j)
            do k=1,ncont_tmp(j)
               if(cont_tmp(k,j).gt.chterm_tmp(1))then
                  cont(k,pos)=cont_tmp(k,j)-chterm_tmp(1)
               else
                  cont(k,pos)=cont_tmp(k,j)-chterm_tmp(1)+chterm_tmp(2)
               endif
            enddo
            !write(*,10)'original:',j,(cont_tmp(k,j),k=1,ncont_tmp(j))
            !write(*,10)'swapped: ',pos,(cont(k,pos),k=1,ncont_tmp(j))
         enddo
      enddo
 10   format(A10,I4,10I3)

      ich=0
      chterm_tmp(0)=0  !boundary condition
      do i=nchain,1,-1
         ich=ich+1
         chname(ich)=chname_tmp(i)
         chterm(ich)=chterm_tmp(i)-chterm_tmp(i-1)+chterm(ich-1)
      enddo

      !get a new mapping from index to chain id
      do i=1,nchain
         do j=chterm(i-1)+1,chterm(i)
            id_chain(j)=i
         enddo
      enddo

!      write(*,'A30')'After exchanging chains...'
!      do i=1,nchain
!         write(*,'(A8,I3)')'Chain ',i
!         do j=chterm(i-1)+1,chterm(i)
!            print *,j,(cont(k,j),k=1,ncont(j))
!         enddo
!      enddo

      return
      end



cccc========================================================cccc
cccc  Swap two PDB structures                               cccc
cccc========================================================cccc
      subroutine swap_structs(coor,nchain1,chterm1,chname1,seq1,pdbid1,
     &     respdbnm1,respdbid1,sec1,id_chain1,ncont1,cont1,lstruct1,
     &     nchain2,chterm2,chname2,seq2,pdbid2,
     &     respdbnm2,respdbid2,sec2,id_chain2,ncont2,cont2,lstruct2)

      include 'pars.h'

      character*5 respdbid1(maxr),respdbid2(maxr), respdbid_tmp(maxr)
      character*3 respdbnm1(maxr),respdbnm2(maxr), respdbnm_tmp(maxr)
      character   chname1(maxk),chname2(maxk), chname_tmp(maxk)
      character   seq1(maxr),seq2(maxr), seq_tmp(maxr)

      real coor(3,maxr,0:1),coor_tmp(3,maxr)

      integer nchain1,nchain2,nchain_tmp
      integer chterm1(0:maxk),chterm2(0:maxk),chterm_tmp(0:maxk)
      integer id_chain1(maxr),id_chain2(maxr)
      integer sec1(maxr),sec2(maxr), sec_tmp(maxr)
      integer ncont1(maxr),ncont2(maxr),ncont_tmp(maxr)
      integer cont1(maxc,maxr),cont2(maxc,maxr), cont_tmp(maxc,maxr)
      integer pdbid1,pdbid2


      !store all variables to temp places
      do i=1,chterm1(nchain1)
         respdbid_tmp(i) = respdbid1(i)
         respdbnm_tmp(i) = respdbnm1(i)
         seq_tmp(i) = seq1(i)

         !print "(A4,I4,1X,A1,A6)",'i= ',i,seq1(i),respdbnm1(i)
         do j=1,3
            coor_tmp(j,i) = coor(j,i,pdbid1)
         enddo

         sec_tmp(i)=sec1(i)
         ncont_tmp(i)=ncont1(i)
         do j=1,ncont1(i)
            cont_tmp(j,i)=cont1(j,i)
         enddo
      enddo

      do i=1,nchain1
         chname_tmp(i)=chname1(i)
         chterm_tmp(i)=chterm1(i)
      enddo

      nchain_tmp=nchain1

      !now swapping, first structure----------------->
      do i=1,chterm2(nchain2)
         respdbid1(i) = respdbid2(i)
         respdbnm1(i) = respdbnm2(i)
         seq1(i) = seq2(i)

         do j=1,3
            coor(j,i,pdbid1) = coor(j,i,pdbid2)
         enddo

         sec1(i)=sec2(i)
         ncont1(i)=ncont2(i)
         do j=1,ncont2(i)
            cont1(j,i)=cont2(j,i)
         enddo
      enddo

      lstruct1=chterm2(nchain2)
      do i=1,nchain2
         chname1(i)=chname2(i)
         chterm1(i)=chterm2(i)
      enddo

      nchain1=nchain2
      do i=1,nchain1
         do j=chterm1(i-1)+1,chterm1(i)
            id_chain1(j)=i
         enddo
      enddo

      !second structure-------------------->
      do i=1,chterm_tmp(nchain_tmp)
         respdbid2(i) = respdbid_tmp(i)
         respdbnm2(i) = respdbnm_tmp(i)
         seq2(i) = seq_tmp(i)

         do j=1,3
            coor(j,i,pdbid2) = coor_tmp(j,i)
         enddo

         sec2(i)=sec_tmp(i)
         ncont2(i)=ncont_tmp(i)
         do j=1,ncont_tmp(i)
            cont2(j,i)=cont_tmp(j,i)
         enddo
      enddo

      lstruct2=chterm_tmp(nchain_tmp)
      do i=1,nchain_tmp
         chname2(i)=chname_tmp(i)
         chterm2(i)=chterm_tmp(i)
      enddo

      nchain2=nchain_tmp
      do i=1,nchain2
         do j=chterm2(i-1)+1,chterm2(i)
            id_chain2(j)=i
         enddo
      enddo


      return
      end


ccccc================================================ccccc
ccccc print the mapping array for debugging purpose  ccccc
ccccc================================================ccccc
      subroutine print_invmap(invmap,len)
      implicit none
      include 'pars.h'

      integer invmap(maxr),len
      integer i,j,prev_i,prev_j

      print *,'printing invmap...'

      prev_i=-1
      prev_j=-1
      do j=1,len
         i=invmap(j)
         if(i.gt.0)then
            if(i-prev_i.ne.1.or.j-prev_j.ne.1)then
               write(*,'(A)') '--'
            endif
            write(*,'(2I6)') i, j
         endif
         prev_i=i
         prev_j=j
      enddo
      write(*,'(/)')

      return
      end



ccccc================================================ccccc
ccccc   reverse invmap 
ccccc================================================ccccc
      subroutine rev_invmap(invmap_r,invmap,nseq1,nseq2)
      implicit none
      include 'pars.h'

      integer invmap(maxr),invmap_r(maxr)
      integer i,j,nseq1,nseq2

      do i=1,nseq2
         invmap(i)=-1
      enddo

      do i=1,nseq1
         j=invmap_r(i)
         if(j.gt.0) invmap(j)=i            !seq2 to seq1
      enddo

      return
      end



ccccc================================================ccccc
ccccc get the length of aligned residues from        ccccc
ccccc invmap                                         ccccc
ccccc================================================ccccc
      subroutine get_alnlen(invmap,nseq2,len)
      implicit none
      include 'pars.h'

      integer nseq2,len,j
      integer invmap(maxr)

      len=0
      do 10 j=1,nseq2
         if(invmap(j).gt.0)then
            len=len+1
         endif
 10   continue

      return
      end

ccccc================================================ccccc
ccccc generate the alignment in sequences for output ccccc
ccccc================================================ccccc
      subroutine get_alnseq(alnmap1,alnmap2,laln,
     &     seq1,seq2,nchain1,nchain2,cterm1,cterm2,
     &     alnx1,alny1,alnz1,alnx2,alny2,alnz2,d_cutoff,
     &     alnseq1,alnseq2,markseq,lseq)

      include 'pars.h'

      character seq1(maxr),seq2(maxr),aa
      character fseq1(maxr),fseq2(maxr)
      character alnseq1(maxr),alnseq2(maxr)
      character markseq(maxr)

      real alnx1(maxr),alny1(maxr),alnz1(maxr)
      real alnx2(maxr),alny2(maxr),alnz2(maxr)
      real d_cutoff  ! distance cutoff for marks

      integer alnmap1(maxr),alnmap2(maxr)
      integer cterm1(0:maxk),cterm2(0:maxk)
      integer laln,lseq
      integer pos    !indexes on final aligned sequence


      lseq=cterm1(nchain1)+cterm2(nchain2)-laln

      do i=1,lseq
         alnseq1(i) = '-'
         alnseq2(i) = '-'
         markseq(i) = ' '
      enddo

      call formatSeq(nchain1,cterm1,seq1,fseq1)
      call formatSeq(nchain2,cterm2,seq2,fseq2)

      pos=0     
      i1_prev=0  !position in unaligned sequence
      i2_prev=0

      alnmap1(laln+1)=cterm1(nchain1)+1  ! artificial alignment for handling sequence ends properly
      alnmap2(laln+1)=cterm2(nchain2)+1
      do 200 k=1,laln+1
         i1=alnmap1(k)
         i2=alnmap2(k)

         do j=i1_prev+1,i1-1
            pos=pos+1
            alnseq1(pos)=fseq1(j)
         enddo

         do j=i2_prev+1,i2-1
            pos=pos+1
            alnseq2(pos)=fseq2(j)
         enddo

         if(k.eq.laln+1)goto 200  !ignore the last artificial alignment

         pos=pos+1
         alnseq1(pos)=fseq1(i1)
         alnseq2(pos)=fseq2(i2)

         dist=sqrt( (alnx1(k)-alnx2(k))**2
     &        +     (alny1(k)-alny2(k))**2
     &        +     (alnz1(k)-alnz2(k))**2 )

         if(dist.le.d_cutoff)then
            markseq(pos)=':'
         else
            markseq(pos)='.'
         endif

         i1_prev=i1
         i2_prev=i2
 200  continue

      return
      end


ccccc==================================================ccccc
ccccc generate sequential order independent alignment  ccccc
ccccc==================================================ccccc
      subroutine get_matchlst(nseq1,nseq2,mf1,mf2,xtmf1,xtmf2,
     &     ytmf1,ytmf2,ztmf1,ztmf2,respdbid1,respdbid2,
     &     seq1,seq2,id_chain1,id_chain2,ncont1,ncont2,
     &     chname1,chname2,naln,matchlst,iter)

      implicit none
      include 'pars.h'
      common/d0/d0,anseq

      character*5 respdbid1(maxr),respdbid2(maxr)
      character   seq1(maxr),seq2(maxr)
      character*3 resnm1,resnm2
      character   chname1(maxk), chname2(maxk)
      character*80 matchlst(maxr)
      character*3  note

      integer mf1(maxr),mf2(maxr)
      real xtmf1(maxr),ytmf1(maxr),ztmf1(maxr)
      real xtmf2(maxr),ytmf2(maxr),ztmf2(maxr)
      real d0,anseq

      integer id_chain1(maxr),id_chain2(maxr)  !original index to chain ID
      integer ncont1(maxr),ncont2(maxr)      !number of contacts
      integer nseq1, nseq2, naln
      integer is_ali(maxr, maxr)   !alignment map
      integer col,ncol
      real    fcol,tms,tms1,ss,dis
      integer iter,i,j,k,i1,i2,k1,k2

      call init_ali_map( is_ali, nseq1, nseq2 )
      do k=1,naln
         i=mf1(k)
         j=mf2(k)
         is_ali(i,j) = 1
      enddo

      tms=0
      tms1=0
      do i=1,naln

         i1=mf1(i)
         i2=mf2(i)
         note='   '

         dis=sqrt((xtmf1(i)-xtmf2(i))**2+
     &        (ytmf1(i)-ytmf2(i))**2+(ztmf1(i)-ztmf2(i))**2)
         if(dis.lt.5) note(2:2) = ':'

         call get_fcol_4s_old(i1,i2,is_ali,fcol,col)
         call AA1toAA3(resnm1, seq1(i1))
         call AA1toAA3(resnm2, seq2(i2))

         if(seq1(i1).eq.seq2(i2).and.seq1(i1).ne.'X')then
            note(3:3) = '*'
         endif

         k1=id_chain1(i1)
         k2=id_chain2(i2)
         if(iter.lt.5) then !proteins not swapped
            write(matchlst(i),20)i,chname1(k1),
     &           respdbid1(i1),resnm1,chname2(k2),
     &           respdbid2(i2),resnm2,dis,col,
     &           ncont1(i1),ncont2(i2),note
         else  ! proteins swapped
            write(matchlst(i),20)i,chname2(k2),
     &           respdbid2(i2),resnm2,chname1(k1),
     &           respdbid1(i1),resnm1,dis,col,
     &           ncont1(i1),ncont2(i2),note
         endif
         ss=1/(1+(dis/d0)**2)
         tms=tms+ss
         tms1=tms1+fcol*ss
         !print *,'fcol=',fcol,' ss=',ss

      enddo
      tms=tms/anseq
      tms1=tms1/anseq
!      write(*,'(a,f5.3,a,f5.3)')'TMS = ',tms,
!     &     ', TMS (cont) = ',tms1

 20   format(i5,2x,2(a3,1x,a5,2x,a3,4x),f6.3,i5,i4,i4,a4)

      return
      end


ccccc=========================================================ccccc
ccccc Convert one-letter amino acid code to three-letter code ccccc
ccccc=========================================================ccccc
      subroutine AA1toAA3(resnm3,resnm1)
      implicit none

      character*3 aa(-1:19),resnm3
      data aa/ 'UNK','GLY','ALA','SER','CYS',
     &     'VAL','THR','ILE','PRO','MET',
     &     'ASP','ASN','LEU','LYS','GLU',
     &     'GLN','ARG','HIS','PHE','TYR',
     &     'TRP'/

      character*1 slc(-1:19),resnm1
      data slc/'X','G','A','S','C',
     &     'V','T','I','P','M',
     &     'D','N','L','K','E',
     &     'Q','R','H','F','Y',
     &     'W'/

      integer j

      resnm3='UNK'
      do j=-1,19
         if(resnm1.eq.slc(j))then
            resnm3=aa(j)
            return
         endif
      enddo

      return
      end


ccccc=========================================================ccccc
ccccc Format sequence: odd chains are lower case, even upper  ccccc
ccccc=========================================================ccccc

      subroutine formatSeq(nchain,cterm,seq,fseq)
      include 'pars.h'

      character*1 seq(maxr),fseq(maxr),aa
      integer cterm(0:maxk)
      integer nchain,ich,i,j

      ! make a copy of the original sequence
      do i=1,maxr
         fseq(i)=seq(i)
      enddo

      ! first make the sequence all upper case
      do i=1,nchain
         do j=cterm(i-1)+1, cterm(i)
            ich=ichar(seq(j))
            if(ich.gt.ichar('Z'))then
               aa =char(ich-32)
               fseq(j)=aa
            endif
         enddo
      enddo

      ! then make sequence of even chains lower case
      do i=1,nchain
         do j=cterm(i-1)+1, cterm(i)
            if(mod(i,2).eq.0)then
               ich=ichar(seq(j))
               aa =char(ich+32)
               fseq(j)=aa
            endif
         enddo
      enddo

      return
      end
ccccc==========================================================cccc


      subroutine getBasename(fullname,basename)
      implicit none
      character fullname*(*),basename*(*)
      integer i,j,l,n_sta,n_end

      l = len(fullname)
      n_sta=0
      n_end=l
      do i=1,l
         if(fullname(i:i).eq.'/')then
            n_sta=i
         else if(fullname(i:i).ne.' ')then
            n_end=i
         endif
      enddo

      l = n_end - n_sta
      do i=1,l
         j=n_sta+i
         basename(i:i)=fullname(j:j)
      enddo

      ! maximum length of basename is 20 characters
      do i=l+1,20
         basename(i:i)=' '
      enddo

      return
      end



ccccc=========================================================ccccc
ccccc Get the transformation matrix for final output          ccccc
ccccc=========================================================ccccc

      subroutine get_TransMatrix(aln2ind,coor,xtrans,ytrans,ztrans,
     &     lenaln,trans,rot)
      include 'pars.h'

      integer aln2ind(maxr),lenaln
      integer L,i,k,ier

      real    coor(3,maxr,0:1),rmsd
      real    xtrans(maxr),ytrans(maxr),ztrans(maxr)
      real*8  trans(3),rot(3,3),rms
      real*8  r_1(3,maxr),r_2(3,maxr),w(maxr)

      data w /maxr*1.0/


      if(lenaln.le.3) return

      do i=1,lenaln
         k=aln2ind(i)
         r_1(1,i)=coor(1,k,0)
         r_1(2,i)=coor(2,k,0)
         r_1(3,i)=coor(3,k,0)
         r_2(1,i)=xtrans(i)
         r_2(2,i)=ytrans(i)
         r_2(3,i)=ztrans(i)
      enddo

      call u3b(w,r_1,r_2,lenaln,1,rms,rot,trans,ier) !u rotate r_1 to r_2
      rmsd=dsqrt(rms/lenaln)

      return
      end

ccccc===================================================ccccc
ccccc  Get inverse transformation matrix.               ccccc
ccccc  The subroutine exploits the property that the    ccccc
ccccc  rotation matrix from RMSD fit is orthogonal, ie, ccccc
ccccc  inverse U = transpose U                          ccccc
ccccc===================================================ccccc
      subroutine inv_TransMatrix( t, u )

      real*8 u(3,3), v(3,3) 
      real*8 t(3), s(3)   

      ! orthogonal matrix, U(inverse) = U(transpose)
      do i=1,3
         do j=1,3
            v(i,j) = u(j,i)
         enddo
      enddo

      ! translation = - U^-1 * t, not -t !
      do i=1,3
         s(i)=0
         do j=1,3
            s(i) = s(i) - v(i,j)*t(j)
         enddo
      enddo

      ! overwrite old matrix
      do i=1,3
         do j=1,3
            u(i,j) = v(i,j)
         enddo
         t(i) = s(i)
      enddo

      return     
      end 

ccccc===================================================ccccc
ccccc  To calculate the radius of gyration              ccccc
ccccc===================================================ccccc
      subroutine calc_rgyr( rg,coor,molid,nseq )
      implicit none
      include 'pars.h'

      real coor(3,maxr,0:1)
      real*8 mx,my,mz,rg

      integer molid,nseq, i

      mx=0
      my=0
      mz=0

      do i = 1, nseq
         mx = mx + coor(1,i,molid)
         my = my + coor(2,i,molid)
         mz = mz + coor(3,i,molid)
      enddo

      mx=mx/nseq
      my=my/nseq
      mz=mz/nseq

      rg=0
      do i = 1,nseq
         rg = rg + (coor(1,i,molid) - mx)**2
         rg = rg + (coor(2,i,molid) - my)**2
         rg = rg + (coor(3,i,molid) - mz)**2
      enddo

      rg=dsqrt(rg/nseq)

      return
      end
      


ccccc===================================================ccccc
ccccc  P-value calculation for IS/TM-score.             ccccc
ccccc  Parameters were obtained empirically from one    ccccc
ccccc  million random interface alignments              ccccc
ccccc===================================================ccccc
      subroutine calcPvalue(score,score_flag,order_flag,nseq1,nseq2,pvalue,z)
      implicit none

      real*8  loc, scale, logt, logq
      real*8  z, pvalue
      real    score
      integer tlen, qlen, nseq1, nseq2
      integer score_flag, order_flag

      qlen = max(nseq1,nseq2)   ! length of the template 
      tlen = min(nseq1,nseq2)   ! length of the target

      if(tlen.lt.20)tlen=20
      if(qlen.lt.20)qlen=20

      logt = log(dble(tlen))
      logq = log(dble(qlen))

      pvalue=-1
      z=-1

      if(order_flag.eq.1)then
         if(score_flag.eq.1)then  ! IS-score
            if( tlen .lt. 45 )then
               loc   = -0.0242 + 0.0391*logt + 0.0071*logq
               scale =  0.0342 - 0.0058*logt + 0.0013*logq
            else 
               loc   = 0.1635  - 0.0085*logt + 0.0057*logq
               scale = 0.0440  - 0.0076*logt + 0.0006*logq
            endif
         else  ! TM-score
            if( tlen.lt.45 )then
               loc   = -0.0690  + 0.0460*logt + 0.0260*logq
               scale =  0.0355  - 0.0050*logt + 0.0024*logq
            else
               loc   =  0.2092 - 0.0368*logt + 0.0344*logq
               scale =  0.0416 - 0.0042*logt + 0.0007*logq
            endif
         endif
      else
         if(score_flag.eq.1)then  ! IS-score
            if( tlen .lt. 55 )then
               loc   = 0.1776 - 0.0038*logt + 0.0113*logq
               scale = 0.0397 - 0.0031*logt - 0.0009*logq
            else 
               loc   = 0.2017 - 0.0160*logt + 0.0163*logq
               scale = 0.0432 - 0.0032*logt - 0.0013*logq
            endif
         else
            if( tlen .lt. 55 )then
               loc   = 0.1986 - 0.0168*logt + 0.0389*logq
               scale = 0.0083 + 0.0097*logt - 0.0011*logq
            else 
               loc   = 0.2398 - 0.0493*logt + 0.0577*logq
               scale = 0.0480 + 0.0023*logt - 0.0028*logq
            endif
         endif
      endif

      if(scale.lt.0.002) scale=0.002 !prevent insane values

      call calcEVDPV(score,loc,scale,pvalue,z)

      !write(*,'(F8.5,2I6,1x,F6.2,g12.4E3)')score,qlen,tlen,z,pvalue
      !write(*,'(2(A,F8.5))')'loc = ',loc,' scale = ',scale
      return     
      end 



ccccc===================================================ccccc
ccccc  P-value calculation for Gumbel Distribution.     ccccc
ccccc===================================================ccccc
      subroutine calcEVDPV(score,loc,scale,pvalue,z)

      real*8  loc, scale
      real*8  z, pvalue
      real    score

      z = (score - loc) / scale
      if( z .lt. 35 ) then
         pvalue = 1 -  dexp( - dexp( -z ) ) !Gumbel distribution
      else
c        below the double precision limit, use taylor expansion approaximation
         pvalue = dexp( -z )
      endif

      return     
      end 


ccccc===================================================ccccc
ccccc  Print the help message                           ccccc
ccccc===================================================ccccc
      subroutine printHelp()

         write(*,*)
         write(*,'(A)')'Usage: IS-align <options> pdbfile1 pdbfile2 '//
     &        'contactfile1 contactfile2'
         write(*,*)
         write(*,'(A)')'Options:'
         write(*,'(A)')'     -t        Use TM-score as the similarity'//
     &        ' measure. The default is IS-score.'

         write(*,'(A)')'     -s        Enable non-sequential alignment.'//
     &        ' By default the alignment is sequential.'
         write(*,'(A)')'     -q <num>  Speed mode:'//
     &     ' 1 - regular (default), 2 - fast (less accurate).'

         write(*,'(A)')'     -L <num>  Normalize the score with a fixed '//
     &        'length specified by num.'

         write(*,'(A)')'     -a        Normalize the score by the average'//
     &        ' size of two interfaces.'
         write(*,'(A)')'     -b        Normalize the score by the minimum'//
     &        ' size of two interfaces.'
         write(*,'(A)')'     -c        Normalize the score by the maximum'//
     &        ' size of two interfaces.'
         write(*,'(A)')'     -v        0 - no alignment, '//
     &     ' 1 - concise alignment, 2 - detailed alignment.'
         write(*,*)

      return
      end




********************************************************************
*     Sequence order dependent matching
*
*     Use an efficient algorithm for solving a linear sum assignment
*     problem at time complexity of O(n^3).
*
*     Mu Gao
********************************************************************
      SUBROUTINE solvLSAP(nseq,score,invmap,nseq1,nseq2,n,debug)
      implicit none
      include 'pars.h'

      integer invmap(maxr)
      integer nseq1, nseq2, nseq
      integer n
      integer spalte(n)
      integer i,j

      real c(n,n), cost, z
      real eps,sup
      real score(nseq,nseq)

      logical debug
      character*50 filename

      eps = 1E-10
      sup = 1E+9


      ! cost matrix for linear sum assignment problem
      do i = 1, n
        do j = 1, n
           c(j,i) = 2
        enddo
      enddo

      do i = 1, n
        do j = 1, n
           if(i.le.nseq1.and.j.le.nseq2)then
              c(i,j) = 2 - score(i,j)
           endif 
        enddo
      enddo

      if(debug)then
         write(filename,"(A,I2,A)")'cmat',nseq2,'.dat'
         open(unit=55,file=filename,status='unknown')
         !write(55,'(2(a,I6))')'len1=',nseq1,' len2=',nseq2
         do j = 1, n
            write(55,'(I6)',advance="no")j
            do i = 1, n
               write(55,'(f11.8)',advance="no")c(j,i)
            enddo
            write(55,*)
         enddo
         close(55)
      endif

      call lsapr ( n, c, z, spalte, sup, eps )

      do j = 1, nseq2
         i=spalte(j)
         if(i.le.nseq1)then
            invmap(j) = i
         else
            invmap(j) = -1
         endif
      end do

      if(debug)write(*,'(a,f11.7)')'optimal cost z=',z

      return
      END
*****************  End of solvLSAP subroutein **************************




      subroutine lsapr ( n, cost, z, spalte, sup, eps )

c*********************************************************************72
c
cc LSAPR solves the linear sum assignment problem with real data.
c
c  Modified:
c
c    22 November 2007
c
c  Author:
c
c    Ulrich Derigs
c
c  Modified by Mu Gao
c
c  Reference:
c
c    Rainer Burkard, Ulrich Derigs,
c    Assignment and Matching Problems: Solution Methods with Fortran-Programs,
c    Springer, 1980,
c    ISBN: 3-540-10267-1,
c    LC: QA402.5 B86.
c
c  Parameters:
c
c    Input, integer N, the dimension of the cost matrix.
c
c    Input, real SUP, a large machine number.
c
c    Input, real C(N,N), the cost matrix.
c
c    Output, real Z, the optimal value.
c
c    Output, integer SPALTE(N), the optimal assignment.
c
c    Input, real EPS, machine accuracy.
c

      implicit none

      integer n

      real cost(n,n)
      real cc
      real d    !Mu: this should be real!
      real dminus(n)
      real dplus(n)
      real eps,sup
      real ui
      real vgl
      real vj
      real ys(n)
      real yt(n)
      real z

      integer i,j,u,w
      integer ind,index,j0
      integer zeile(n)
      integer vor(n)
      integer spalte(n)

      logical label(n)

c
c  Construct an initial partial assignment.
c
      do i = 1, n
        zeile(i) = 0
        spalte(i) = 0
        vor(i) = 0
        ys(i) = 0.0E+00
        yt(i) = 0.0E+00
      end do

      do i = 1, n
        do j = 1, n
          cc = cost(j,i)

          if ( j .ne. 1 ) then
            if ( ( cc - ui ) .ge. eps ) then
              go to 3
            end if
          end if

          ui = cc
          j0 = j

3         continue
        end do

        ys(i) = ui

        if ( zeile(j0) .eq. 0 ) then
          zeile(j0) = i
          spalte(i) = j0
        end if

      end do

      do j = 1, n
        yt(j) = 0
        if ( zeile(j) .eq. 0 ) then
          yt(j) = sup
        end if
      end do


      do i = 1, n
        ui = ys(i)

        do j = 1, n
          vj = yt(j)
          if ( eps .lt. vj ) then

            cc = cost(j,i) - ui
            if ( cc + eps .lt. vj ) then
              yt(j) = cc
              vor(j) = i
            end if

          end if
        end do
      end do

      do j = 1, n
        i = vor(j)
        if ( i .ne. 0 ) then
          if ( spalte(i) .eq. 0 ) then
            spalte(i) = j
            zeile(j) = i
          end if
        end if
      end do

      do i = 1, n
        if ( spalte(i) .eq. 0 ) then
          ui = ys(i)

          do j = 1, n
            if ( zeile(j) .eq. 0 ) then
              cc = cost(j,i)
              if ( cc - ui - yt(j) .le. -eps ) then
              !if ( ( cc - ui - yt(j) + eps ) .le. 0.0E+00 ) then
                spalte(i) = j
                zeile(j) = i
              end if
            end if
          end do

        end if
      end do
c
c  Construct the optimal solution.
c
      do 1000 u = 1, n

        if ( spalte(u) .gt. 0 ) goto 1000

c
c  Shortest path computation.
c

        do i = 1, n
          vor(i) = u
          label(i) = .false.
          dplus(i) = sup
          dminus(i) = cost(i,u) - ys(u) - yt(i)
        end do

        dplus(u) = 0.0E+00

105     continue

        d = sup
        index = 0

        do i = 1, n
          if ( .not. label(i) ) then
            if ( dminus(i) + eps .lt. d ) then
              d = dminus(i)
              index = i
            end if
          end if
        end do

        if ( index .eq. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'LSAPR - Fatal error!'
          write ( *, '(a)' ) '  No unlabeled node with DMINUS < D.'
          stop
        end if

        if ( zeile(index) .le. 0 ) then
          go to 400
        end if

        label(index) = .true.
        w = zeile(index)
        dplus(w) = d

        do i = 1, n
          if ( .not. label(i) ) then
            vgl = d + cost(i,w) - ys(w) - yt(i)
            if ( vgl + eps .lt. dminus(i) ) then
              dminus(i) = vgl
              vor(i) = w
            end if
          end if
        end do

        go to 105
c
c  Augmentation.
c
400     continue

        w = vor(index)
        zeile(index) = w
        ind = spalte(w)
        spalte(w) = index

        if ( w .ne. u ) then
          index = ind
          go to 400
        end if
c
c  Transformation.
c
!500     continue

        do i = 1, n
          if ( dplus(i) .ne. sup ) then
            ys(i) = ys(i) + d - dplus(i)
          end if

          if ( dminus(i) + eps .lt. d ) then
            yt(i) = yt(i) + dminus(i) - d
          end if
        end do

1000  continue
c
c  Computation of the optimal value.
c
      z = 0.0E+00
      do i = 1, n
        j = spalte(i)
        z = z + cost(j,i)
      end do

      return
      end
*****************  End of LSAPR subroutein **************************




cccccc==================================================
ccccc extract a submatrix
ccccc 
ccccc===================================================
      subroutine extract_score(nseq,score,nseq1,nseq2,
     & i1,j1,i2,j2,new_score)
      implicit none
      include 'pars.h'

      real    new_score(nseq,nseq),score(nseq,nseq)
      integer i,j
      integer nseq,nseq1,nseq2,i1,j1,i2,j2

      if(i1.lt.1.or.i1.gt.nseq1)then
         write(*,*)'Error: out of the scoring matrix boundary!'
         stop
      endif

      if(i2.lt.1.or.i2.gt.nseq1)then
         write(*,*)'Error: out of the scoring matrix boundary!'
         stop
      endif
      
      if(j1.lt.1.or.j1.gt.nseq2)then
         write(*,*)'Error: out of the scoring matrix boundary!'
         stop
      endif

      if(j2.lt.1.or.j2.gt.nseq2)then
         write(*,*)'Error: out of the scoring matrix boundary!'
         stop
      endif

      if(i1.ge.i2)then
         write(*,*)'Error: out of the scoring matrix boundary!'
         stop
      endif

      if(j1.ge.j2)then
         write(*,*)'Error: out of the scoring matrix boundary!'
         stop
      endif

      do i=i1,i2
         do j=j1,j2
            new_score(i-i1+1,j-j1+1)=score(i,j)
         enddo
      enddo

      return
      end


cccccc==================================================
ccccc combine two sub invmaps
ccccc===================================================
      subroutine comb_invmap(new_invmap,invmap1,invmap2,
     &     i1,j1,j2,nseq2)
      implicit none
      include 'pars.h'

      integer new_invmap(maxr),invmap1(maxr),invmap2(maxr)
      integer i1,j1,j2,nseq2
      integer i,j

      if(j1+j2.ne.nseq2)then
         write(*,*)'Error: cannot combine two invmaps'
      endif

      do j=1,nseq2
         new_invmap(j)=-1
      enddo

      do j=1,j1
         i=invmap1(j)
         if(i.gt.0) new_invmap(j)=i
      enddo

      do j=1,j2
         i=invmap2(j)
         if(i.gt.0) new_invmap(j+j1)=i+i1
      enddo

      return
      end


cccccc==================================================
ccccc initialize a matrix
ccccc===================================================


      subroutine init_ali_map( is_ali, nseq1, nseq2 )
      include 'pars.h'

      integer nseq1,nseq2
      integer is_ali(maxr, maxr)   !alignment map
      integer i,j

ccc   is_ali(i,j) indicate whether i, j are in alignment dis<d_col
      do i=1,nseq1
         do j=1,nseq2
            is_ali(i,j) = -1
         enddo
      enddo

      return
      end


ccccc==========================================================cccc

      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ', 
     &  'May      ', 'June     ', 'July     ', 'August   ', 
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *, 
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) 
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end



cccc================================================================
