%This script seeks to check for duplicates in Interbase parameter files
%from the G4 Analysis. It takes sets of 8 values and compares them to every
%set of other 8 values.
%NOTE: This algorithm could be more efficient if I only compared lines of
%data with other lines of data. Perhaps, instead of doing all 8 values in a
%column, I could've chosen all 6 values in a row, and just checked for
%every other line that contains that same information since the repeat
%would be found every other 8 files.

IB = readmatrix("3QXR.csv")
Iterations = 0;
for i = 0:998
    st = i*8+1;
    nd = st+7;
    A = IB(st:nd,1);
    for j = 0:998-i
        st2 = nd+1+j*8;
        nd2 = st2+7;
        B = IB(st2:nd2,1);
        if A == B
            fprintf("Duplicate")
            i
            j
        else
            Iterations = Iterations + 1;
        end
    end
end
Iterations