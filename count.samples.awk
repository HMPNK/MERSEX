function logfact(n,   i, sum) {
    if (n <= 1) return 0
    sum = 0
    for (i = 2; i <= n; i++) sum += log(i)
    return sum
}

function fisher_log(a, b, c, d,   n, log_p, p) {
    n = a + b + c + d

    log_p =  logfact(a + b) + logfact(c + d) + logfact(a + c) + logfact(b + d) \
           - (logfact(a) + logfact(b) + logfact(c) + logfact(d) + logfact(n))

    p = exp(log_p)
    return p
}


BEGIN{
if(cmin1==""){cmin1=2};
if(cmin2==""){cmin2=1};
}

{
n=split($0,ar,"\t");
f=0;f2=0;m=0;m2=0;
for(x=2;x<=fcount+1;x++)        {
                        if(ar[x]>=cmin1){f++} else if(ar[x]<=cmin2){f2++}
                        };

for(x=fcount+2;x<=n;x++)        {
                        if(ar[x]>=cmin1){m++} else if(ar[x]<=cmin2){m2++}
                        };

    a = f
    b = f2
    c = m
    d = m2

    p = fisher_log(a, b, c, d)

print p"\t"a"\t"b"\t"c"\t"d"\t"$0;
}
