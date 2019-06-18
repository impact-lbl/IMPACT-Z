
def conciseReadInput(inputFileName):
    try:
        fin = open(inputFileName,'r')
        dataList  = fin.readlines()
        fin.close()
    except:
        print(( "  ERROR! Can't open file '" + inputFileName + "'"))
        return False

    i=0
    while i < len(dataList):
        if dataList[i].lstrip()=='' or dataList[i].lstrip().startswith('!'):
            del dataList[i]
            i=i-1
        else:
            index = dataList[i].lstrip().find('!')
            if index==-1:
                dataList[i]=dataList[i].strip()+'\n'
            else:
                dataList[i]=dataList[i].lstrip()[:index].rstrip()+'\n'
        i=i+1
    dataList  = [line.split() for line in dataList ]
    
    for i in range(0,len(dataList)):
        for j in range(0,len(dataList[i])):
            dataList[i][j] = DtoE(dataList[i][j])
    return dataList


def conciseWriteInput(dataList,outputFileName):
    try:
        fin = open(outputFileName,'w')
    except:
        print(( "  ERROR! Can't open file '" + outputFileName + "'"))
        return False

    for line in dataList:
        for elem in line:
            ImpactInput.writelines(elem+' ')
        ImpactInput.writelines('\n')
    return


def DtoE(word):
    if 'D' in word or 'd' in word: 
        try:
            temp = float(word.replace('D','E',1).replace('d','e',1))
            return str(temp)
        except:
            return word
    else:
        return word