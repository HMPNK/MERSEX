{
n=split($0,a,"\t");

for(x=2;x<=fcount+1;x++)        {
                        if(a[x]>=cmin){f++}
                        };

for(x=fcount+2;x<=n;x++)        {
                        if(a[x]>=cmin){m++}
                        };

print f"\t"m"\t"$0;
f=0;
m=0;
}
