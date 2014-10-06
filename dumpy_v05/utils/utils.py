def str_suffix(n,length=6):
    fewzero=''
    for i in range(length-len(str(n))):
        fewzero=fewzero+'0'
    return fewzero+str(n)
