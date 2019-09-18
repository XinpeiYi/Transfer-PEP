# Input transfer fdr case:
```
% CASE = 1 : fG0 = f0, fG1 = f1
% CASE = 2 : fG0 != f0, fG1 = f1
% CASE = 3 : fG0 = f0, fG1 != f1
% CASE = 4 : fG0 != f0, fG1 != f1
CASE =3;
fdrthres = 0.01; % FDR threshold
```
# Read Mascot Result
```
pathin = '/Users/yixinpei/PosdocResearch/Transfer fdr/code/testData/fG0!=f0_fG1!=f1';
[result] = ReadDatResultFolder(pathin);
```
# Judge Group
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
# Compute Target-Decoy FDR
```
[FDR,Iid,Threshold,FinalFDR,P] = ComputeFDR(DecoyType.total,GroupType.total,scores.total,numrst.total,I.total,fdrthres);
```
# Global fdr
```
ems.global = 0.1;
ppi0 = 0.3;
ppi1 = 0.7;
e = 0.1:0.1:10;
data_global = zeros(1,length(e));
g = 1000000;
for i = 1:length(e)
    i
    hh1.global = e(i);
    [p.global,pi0.global,pi1.global,h0.global,h1.global,f0.global,f1.global] = SemiParametricFitting(scores.decoy,scores.target,ems.global,ppi0,ppi1,hh1.global);
    [F0.global,F1.global] = ComputeF(scores.target,scores.decoy,h0.global,h1.global,p.global);       
    
    FDRfdr.global = (pi0.global*(1-F0.global))./(pi0.global*(1-F0.global)+pi1.global*(1-F1.global));
    
    f = sum((FDR.GF(DecoyType.total==1)-FDRfdr.global).^2);
    
    if f>g
        i-1
        break;
    end
    data_global(i) = f;
    g = f;
end
hh1.global = e(i-1);
ppi0 = 0.3;
ppi1 = 0.7;

[p.global,pi0.global,pi1.global,h0.global,h1.global,f0.global,f1.global] = SemiParametricFitting(scores.decoy,scores.target,ems.global,ppi0,ppi1,hh1.global);
[F0.global,F1.global] = ComputeF(scores.target,scores.decoy,h0.global,h1.global,p.global); 
```
# Separate fdr
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
# Transfer fdr
```
if CASE == 1
    maxiter = 10000;
    f0.trans = f0.global(GroupType.target==0);
    f1.trans = f1.global(GroupType.target==0);
    F0.trans = F0.global(GroupType.target==0);
    F1.trans = F1.global(GroupType.target==0);
    [pi1.trans,pi0.trans] = EM_Change_Pi(f0.trans,f1.trans,scores.target(GroupType.target==0),maxiter);
    [fdrfdr,FDRfdr,Iidfdr,Thresholdfdr,FinalFDRfdr] = ComputelocalFDR_CASE1(pi0,pi1,f0,f1,F0,F1,I.target,fdrthres,GroupType.target,scores.target,numrst.target);
    DrawingThreeFitting(scores.target,GroupType.target,scores.decoy,GroupType.decoy,h0,h1,pi0,pi1,p,f0,f1);
else
    if CASE == 2
        A = -(pi0.global*(P(1))*(1-F0.global(GroupType.target==0))+pi0.global*(P(1)*scores.target(GroupType.target==0)+P(2)).*(-f0.global(GroupType.target==0)));
        pk = sum(GroupType.target==0)/length(scores.target);
        maxiter = 10000;
        II = find(GroupType.target==0); 
        f1.trans = f1.global(II);
        [pi1.trans,pi0.trans] = EM_Difdis_pi(f1.trans,A,pk,scores.target,maxiter);
        f0.trans = A./(pk*pi0.trans);
        scores_group = scores.target(II);
        F1.trans = F1.global(II);
        F0.trans = 1-((P(1)*scores_group+P(2))*pi0.global.*(1-F0.global(II)))./(pi0.trans*pk);
        [fdrfdr,FDRfdr,Iidfdr,Thresholdfdr,FinalFDRfdr] = ComputelocalFDR_CASE2(pi0,pi1,f0,f1,F0,F1,I.target,fdrthres,GroupType.target,scores.target,numrst.target);
        DrawingThreeFitting(scores.target,GroupType.target,scores.decoy,GroupType.decoy,h0,h1,pi0,pi1,p,f0,f1);
    else
        if CASE == 3
            f0.trans = f0.global(GroupType.target==0);
            F0.trans = F0.global(GroupType.target==0);
            
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

            A = -(pi0.global*(P(1))*(1-F0.global(GroupType.target==0))+pi0.global*(P(1)*Targetscores_group+P(2)).*(-f0.global(GroupType.target==0)));
            pk = sum(GroupType.target==0)/length(scores.target);

            for i = 1:length(c)
                i
                hh1.trans = c(i);
                [p.trans,pi0.trans,pi1.trans,h1.trans,f0.trans,f1.trans] = SemiParametricFitting_trans_CASE4(Targetscores_group,ems.trans,ppi0,ppi1,hh1.trans,A,pk);
                [F1.trans,pi0A_1_F0A] = ComputeF_trans_CASE4(Targetscores_group,h1.trans,p.trans,P,pi0.global,F0.global(GroupType.target==0),pk);
                
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
            [F1.trans,pi0A_1_F0A] = ComputeF_trans_CASE4(Targetscores_group,h1.trans,p.trans,P,pi0.global,F0.global(GroupType.target==0),pk);
            [fdrfdr,FDRfdr,Iidfdr,Thresholdfdr,FinalFDRfdr] = ComputelocalFDR_CASE4(pi0,pi1,f0,f1,F0,F1,pi0A_1_F0A,I.target,fdrthres,GroupType.target,scores.target,numrst.target);
            DrawingThreeFitting(scores.target,GroupType.target,scores.decoy,GroupType.decoy,h0,h1,pi0,pi1,p,f0,f1);
        end
    end
end
```
# Drawing the consistency between two FDR estimated methods
```
DrawingFDRfdr(FDR,FDRfdr,DecoyType.total,GroupType.total);
```
