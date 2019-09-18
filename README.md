# Grouped local false discovery rate estimation for peptide identification

Transfer fdr is the first grouped fdr estimation algorithm, which is proposed for quality control of small groups of peptide identifications.
Transfer fdr derives the group null distribution through its empirical relationship with the combined null distribution, and estimates the group alternative distribution, as well as the null proportion, using an iterative semi-parametric method.
Validated on both simulated data and real proteomic data, transfer fdr showed remarkably higher accuracy than the direct combined and separate fdr estimation methods.

Tranfer fdr is implemented in the matlab programming language. The input is the Mascot(\*.dat) identifications. The output is three fdr estimation results: combined fdr, separate fdr and transfer fdr. Here is the main function description.

## Select algorithm case and specify the FDR threshold:
Four different cases are corresponding to four conditions of transfer fdr algorithm, respectively.
Default is case 3.

```
% CASE = 1 : fG0 = f0, fG1 = f1
% CASE = 2 : fG0 != f0, fG1 = f1
% CASE = 3 : fG0 = f0, fG1 != f1
% CASE = 4 : fG0 != f0, fG1 != f1
CASE =3;
fdrthres = 0.01; % FDR threshold
```
## Read Mascot Result
The input is the path of Mascot identifications.
```
pathin = '/Users/yixinpei/PosdocResearch/Transfer fdr/code/testData/fG0!=f0_fG1!=f1';
[result] = ReadDatResultFolder(pathin);
```
## Judge Group
Specify the tags for both the decoy database and the group identifications.
The group identifications can be modification or one certain protein identifications.

```
% GroupType: 0 is Group; 1 is nonGroup.
% DecoyType: 0 is Decoy; 1 is Target.
DecoyTag = 'REVERSE_'; % the tag for the decoy database
TagType = 'Modification'; % can be the modification or one certain protein
GroupTag = 'Methyl'; % the tag for the group

[DecoyType.total,GroupType.total,scores.total,numrst.total,I.total] = JudgeGroup(result,TagType,DecoyTag,GroupTag);

result_sort = result(I.total);
result_target = result_sort(DecoyType.total==1);
result_decoy = result_sort(DecoyType.total==0);
[DecoyType.target,GroupType.target,scores.target,numrst.target,I.target] = JudgeGroup(result_target,TagType,DecoyTag,GroupTag);
[DecoyType.decoy,GroupType.decoy,scores.decoy,numrst.decoy,I.decoy] = JudgeGroup(result_decoy,TagType,DecoyTag,GroupTag);
```
## Compute Target-Decoy FDR

In this part, the three FDR estimation methods: iCombined FDR, iSeparate FDR and iTransfer FDR are calculated.
The output is the three FDR estimations, filtered numbers, filtered threshold and both the slope and intercept for iTransfer FDR. 

```
[FDR,Iid,Threshold,FinalFDR,P] = ComputeFDR(DecoyType.total,GroupType.total,scores.total,numrst.total,I.total,fdrthres);
```
## Combined fdr
Combined fdr is estimated by using iterative semi-parametric method on the whole set of PSMs.

```
ems.combined = 0.1; # Convergence condition
ppi0 = 0.3;
ppi1 = 0.7;
e = 0.1:0.1:10;
data_combined = zeros(1,length(e));
g = 1000000;
for i = 1:length(e)
    i
    hh1.combined = e(i);
    [p.combined,pi0.combined,pi1.combined,h0.combined,h1.combined,f0.combined,f1.combined] = SemiParametricFitting(scores.decoy,scores.target,ems.combined,ppi0,ppi1,hh1.combined);
    [F0.combined,F1.combined] = ComputeF(scores.target,scores.decoy,h0.combined,h1.combined,p.combined);       
    
    FDRfdr.combined = (pi0.combined*(1-F0.combined))./(pi0.combined*(1-F0.combined)+pi1.combined*(1-F1.combined));
    
    f = sum((FDR.GF(DecoyType.total==1)-FDRfdr.combined).^2);
    
    if f>g
        i-1
        break;
    end
    data_combined(i) = f;
    g = f;
end
hh1.combined = e(i-1);
ppi0 = 0.3;
ppi1 = 0.7;

[p.combined,pi0.combined,pi1.combined,h0.combined,h1.combined,f0.combined,f1.combined] = SemiParametricFitting(scores.decoy,scores.target,ems.combined,ppi0,ppi1,hh1.combined);
[F0.combined,F1.combined] = ComputeF(scores.target,scores.decoy,h0.combined,h1.combined,p.combined); 
```
## Separate fdr
Separate fdr is estimated by using iterative semi-parametric method on the PSMs in the group only.

```
ppi0 = 0.3;
ppi1 = 0.7;
ems.separate = 0.05;
Targetscores_group = scores.target(GroupType.target==0);
Decoyscores_group = scores.decoy(GroupType.decoy==0);
c = 0.5:0.1:10;
data_separate = zeros(1,length(c));
d = 1000000;
for i = 1:length(c)
    i
    hh1.separate = c(i);
    [p.separate,pi0.separate,pi1.separate,h0.separate,h1.separate,f0.separate,f1.separate] = SemiParametricFitting(Decoyscores_group,Targetscores_group,ems.separate,ppi0,ppi1,hh1.separate);
    [F0.separate,F1.separate] = ComputeF(Targetscores_group,Decoyscores_group,h0.separate,h1.separate,p.separate);    
    
    FDRfdr.separate = (pi0.separate*(1-F0.separate))./(pi0.separate*(1-F0.separate)+pi1.separate*(1-F1.separate));
    
    b = sum((FDR.SF(DecoyType.total(GroupType.total==0)==1)-FDRfdr.separate).^2);
    
    if b>d
        i-1
        break;
    end
    data_separate(i) = b;
    d = b;
end
hh1.separate = c(i-1);
[p.separate,pi0.separate,pi1.separate,h0.separate,h1.separate,f0.separate,f1.separate] = SemiParametricFitting(Decoyscores_group,Targetscores_group,ems.separate,ppi0,ppi1,hh1.separate);
[F0.separate,F1.separate] = ComputeF(Targetscores_group,Decoyscores_group,h0.separate,h1.separate,p.separate);    
```
## Transfer fdr
Transfer fdr is estimated based on the case selected in the first part.

```
if CASE == 1
    maxiter = 10000;
    f0.trans = f0.combined(GroupType.target==0);
    f1.trans = f1.combined(GroupType.target==0);
    F0.trans = F0.combined(GroupType.target==0);
    F1.trans = F1.combined(GroupType.target==0);
    [pi1.trans,pi0.trans] = EM_Change_Pi(f0.trans,f1.trans,scores.target(GroupType.target==0),maxiter);
    [fdrfdr,FDRfdr,Iidfdr,Thresholdfdr,FinalFDRfdr] = ComputelocalFDR_CASE1(pi0,pi1,f0,f1,F0,F1,I.target,fdrthres,GroupType.target,scores.target,numrst.target);
    DrawingThreeFitting(scores.target,GroupType.target,scores.decoy,GroupType.decoy,h0,h1,pi0,pi1,p,f0,f1);
else
    if CASE == 2
        A = -(pi0.combined*(P(1))*(1-F0.combined(GroupType.target==0))+pi0.combined*(P(1)*scores.target(GroupType.target==0)+P(2)).*(-f0.combined(GroupType.target==0)));
        pk = sum(GroupType.target==0)/length(scores.target);
        maxiter = 10000;
        II = find(GroupType.target==0); 
        f1.trans = f1.combined(II);
        [pi1.trans,pi0.trans] = EM_Difdis_pi(f1.trans,A,pk,scores.target,maxiter);
        f0.trans = A./(pk*pi0.trans);
        scores_group = scores.target(II);
        F1.trans = F1.combined(II);
        F0.trans = 1-((P(1)*scores_group+P(2))*pi0.combined.*(1-F0.combined(II)))./(pi0.trans*pk);
        [fdrfdr,FDRfdr,Iidfdr,Thresholdfdr,FinalFDRfdr] = ComputelocalFDR_CASE2(pi0,pi1,f0,f1,F0,F1,I.target,fdrthres,GroupType.target,scores.target,numrst.target);
        DrawingThreeFitting(scores.target,GroupType.target,scores.decoy,GroupType.decoy,h0,h1,pi0,pi1,p,f0,f1);
    else
        if CASE == 3
            f0.trans = f0.combined(GroupType.target==0);
            F0.trans = F0.combined(GroupType.target==0);
            
            ppi0 = 0.3;
            ppi1 = 0.7;
            ems.trans = 0.05;
            c = 0.5:0.1:10;
            data_trans = zeros(1,length(c));
            d = 1000000000;
            
            for i = 1:length(c)
                i
                hh1.trans = c(i);
                [p.trans,pi0.trans,pi1.trans,h1.trans,f0.trans,f1.trans] = SemiParametricFitting_trans_CASE3(Targetscores_group,ems.trans,ppi0,ppi1,hh1.trans,f0.trans);
                [F1.trans] = ComputeF_trans_CASE3(Targetscores_group,h1.trans,p.trans);
                
                [FDRfdr.trans] = ComputelocalFDR_trans_CASE3(pi0.trans,pi1.trans,F0.trans,F1.trans);
                
                b = sum((FDR.TF(DecoyType.total(GroupType.total==0)==1)-FDRfdr.trans).^2);
                
                if b>d
                   i-1
                   break;
                end
                data_trans(i) = b;
                d = b;
            end
            hh1.trans = c(i-1);
            [p.trans,pi0.trans,pi1.trans,h1.trans,f0.trans,f1.trans] = SemiParametricFitting_trans_CASE3(Targetscores_group,ems.trans,ppi0,ppi1,hh1.trans,f0.trans);
            [F1.trans] = ComputeF_trans_CASE3(Targetscores_group,h1.trans,p.trans);
            [fdrfdr,FDRfdr,Iidfdr,Thresholdfdr,FinalFDRfdr] = ComputelocalFDR_CASE3(pi0,pi1,f0,f1,F0,F1,I.target,fdrthres,GroupType.target,scores.target,numrst.target);
            DrawingThreeFitting(scores.target,GroupType.target,scores.decoy,GroupType.decoy,h0,h1,pi0,pi1,p,f0,f1);
        else
            ppi0 = 0.3;
            ppi1 = 0.7;
            ems.trans = 0.05;
            c = 0.5:0.1:10;
            data_trans = zeros(1,length(c));
            d = 1000000000;

            A = -(pi0.combined*(P(1))*(1-F0.combined(GroupType.target==0))+pi0.combined*(P(1)*Targetscores_group+P(2)).*(-f0.combined(GroupType.target==0)));
            pk = sum(GroupType.target==0)/length(scores.target);

            for i = 1:length(c)
                i
                hh1.trans = c(i);
                [p.trans,pi0.trans,pi1.trans,h1.trans,f0.trans,f1.trans] = SemiParametricFitting_trans_CASE4(Targetscores_group,ems.trans,ppi0,ppi1,hh1.trans,A,pk);
                [F1.trans,pi0A_1_F0A] = ComputeF_trans_CASE4(Targetscores_group,h1.trans,p.trans,P,pi0.combined,F0.combined(GroupType.target==0),pk);
                
                [FDRfdr.trans] = ComputelocalFDR_trans_CASE4(pi1.trans,pi0A_1_F0A,F1.trans);
                
                b = sum((FDR.TF(DecoyType.total(GroupType.total==0)==1)-FDRfdr.trans).^2);
    
                if b>d
                   i-1
                   break;
                end
                data_trans(i) = b;
                d = b;
            end
            hh1.trans = c(i-1);
            [p.trans,pi0.trans,pi1.trans,h1.trans,f0.trans,f1.trans] = SemiParametricFitting_trans_CASE4(Targetscores_group,ems.trans,ppi0,ppi1,hh1.trans,A,pk);
            [F1.trans,pi0A_1_F0A] = ComputeF_trans_CASE4(Targetscores_group,h1.trans,p.trans,P,pi0.combined,F0.combined(GroupType.target==0),pk);
            [fdrfdr,FDRfdr,Iidfdr,Thresholdfdr,FinalFDRfdr] = ComputelocalFDR_CASE4(pi0,pi1,f0,f1,F0,F1,pi0A_1_F0A,I.target,fdrthres,GroupType.target,scores.target,numrst.target);
            DrawingThreeFitting(scores.target,GroupType.target,scores.decoy,GroupType.decoy,h0,h1,pi0,pi1,p,f0,f1);
        end
    end
end
```
## Drawing the consistency between two FDR estimated methods


```
DrawingFDRfdr(FDR,FDRfdr,DecoyType.total,GroupType.total);
```


# Questions and Technical Support
If your have any questions, comments, suggestions or problems, please let us know.

For more informations about transfer fdr algorithm, e.g. test data or references, see our website http://fugroup.amss.ac.cn/software/TransferFdr/TransferFdr.html.
