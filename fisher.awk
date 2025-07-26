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
    # Beispiel: 2x2-Tabelle mit festen Werten
    a = $1
    b = nf-$1
    c = $2
    d = nm-$2

print a"\t"b
print c"\t"d

    p = fisher_log(a, b, c, d)
    print p"\t"$0;
}
