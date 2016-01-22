%-1 is empty and not saved
%0 and 1 are saved
function saveParts(A, prefix)
B = A;
B( abs(B)<0.5) = 0;
B( B~=0 ) = 1;
saveArrInt2d(1-B, strcat(prefix,'0.txt'))
B = A;
B(B<0.5) = 0;
B(B>0.5) = 1;
saveArrInt2d(B, strcat(prefix, '1.txt'));
end