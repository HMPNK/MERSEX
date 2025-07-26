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

{
n=split($0,a,"\t");
f=0;f2=0;m=0;m2=0;
for(x=2;x<=fcount+1;x++)        {
                        if(a[x]>=cmin){f++} else{f2++}
                        };

for(x=fcount+2;x<=n;x++)        {
                        if(a[x]>=cmin){m++} else{m2++}
                        };

    a = f
    b = f2
    c = m2
    d = m

    p = fisher_log(a, b, c, d)

print p"\t"a"\t"b"\t"c"\t"d"\t"$0;
}
