function pR = ChooseMixedShredLaw(Wr, Ws )
%ChooseMixedShredLaw(Wr, Ws ) Determine which shred law will be mixed with
%the remainder law comparing the weight Wr of the remainder law and the required weight of each law.

pR=1;
nl=length(Ws);
while( pR<= nl && Ws(pR) < Wr  )
    pR=pR+1;
end
if(pR>nl)
    header='ShredByAffinity : \n';
    errmsg=sprintf('The weight of the remainder law (%f) exceed the maximal required weight for the shred laws (%f) \n', Wr, max(Ws) );
    errsug= 'Try increasing the number of points or the maximal required weight';
    error('MJP:InvalidShredingArgument', [header errmsg errsug]);
end

end

