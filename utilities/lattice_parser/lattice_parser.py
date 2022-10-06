import sys
import re
from math import *

class lattice_parser:
    
    def __init__(self, fileName, lineName):
          
        self.fileName = fileName          
        self.lineName = lineName.upper()  # in code, upper case convention
        
    def get_brieflines(self):
        '''
        # delete the comments, blank lines, and so on
        '''
        lines = self._readfile()
        
        lines = self._combine_lines(lines)
        lines = self._delete_comment_blank_lines(lines)
        lines = self._delete_redundant_comma(lines)
        lines = self._delete_backslashn(lines)

        return lines    

    def get_trackline(self, lines):
        '''
        get tracked line.
        lines should be lattice section lines.
        '''
        lattice = self._get_lattice(lines) 

        usedline = lattice[self.lineName]['LINE'].split(',')
        # expand the track line
        usedline = self._expandline(usedline, lines) 
        
        trackline = []
        for item in usedline:
            trackline.append( lattice[item] )
            
        return trackline
    
    # sub-funcs
    #======================================================================
    def _readfile(self):
        f = open(self.fileName,'r')
        lines = f.readlines()
        f.close()

        return lines
    
    def _combine_lines(self,lines):
        '''combine lines ended with & to one line'''

        j = 0
        for line in lines:
            #be careful of "\t": ..., & \t \n
            if re.match(r'.*&+ *\t* *\n$',line):
                # delete '& \n'
                lines[j] = re.sub(r'&+.*\n','',lines[j])
    
                # combine these two line
                lines[j+1] = lines[j]+lines[j+1]
                
                # replace next line with comment line
                # to keep len(lines) unchanged, because in the for loop
                lines[j] = '! \n'
            j = j+1
                
        return lines

    def _delete_comment_blank_lines(self, lines):
        '''
        delete comment (start with !) and blank lines,
        white spaces in each line are also deleted;
        delete comments at the end of each line;              
        '''
        # get the index of comment and blank line
        j = 0
        index=[]
        for line in lines:
            # get the index of comment and blank lines
            if re.match(r'^!',line):
                index.append(j)
            
            # blank lines, pay attention to \t tab space    
            elif re.match(r'^\t* *\n',line):
                index.append(j)
                        
            j = j+1
        
        # delete the comment and blank lines
        cnt = 0
        for j in index:
            del lines[j-cnt]
            cnt += 1
          
        # delete all comments at the end of each line
        j = 0
        for line in lines:
            lines[j] = line.split('!')[0]
            j += 1
            
        return lines

    def _delete_redundant_comma(self, lines):
        '''
        delete all redundant comma,
        delete quatation marks in each line.
        '''
        j = 0
        for line in lines:
            # rpn expression expand, change here
            # remove quatation marks in line, ' ' and " "
            #line = line.replace('\'','').replace('\"','')
    
            # remove redundant comma
            tmplist = line.split(',')
            
            # remove all the '' in tmplist
            while '' in tmplist:
                tmplist.remove('')
    
            # replace lines[j] with no-redundant comma string
            lines[j] = ','.join( tmplist)
            j += 1
        
        return lines

    def _delete_backslashn(self, lines):
        '''
        delete all \n in lines
        '''
        j = 0
        for line in lines:            
            # delete ;\n
            lines[j] = line.replace(';\n','')
         
            # in case no semicolon at the end of each line
            # only \n at the end
            line = lines[j]
            lines[j] = line.replace('\n','')
            
            # \t at the end of line
            line = lines[j]
            lines[j] = line.replace('\t','')
            
            j += 1
        return lines
  
    def _get_lattice(self, lines):
        '''
        save lattice lines in a dict, like:
            {'D1': {'NAME': 'D1', 'TYPE': 'DRIFT', 'L': '1.0'},
            'B1': {'NAME': 'B1', 'TYPE': 'BEND', 'ANGLE': '1.0'}
              ...
            'USELINE': {'NAME': 'USELINE', 'TYPE': 'LINE', 'LINE': 'D1,4*BC'}
            }
        '''
        # lines = self.get_brieflines() #only lattice section
        
        lattice = {}
        for line in lines:
            # there are four types line,
            # 1. L1_q1: quad, L=0.2
            # 2. line: line=(L1_q1)
            # 3. Lq=0.2
            # 4. % 0.2 sto L1=q 
            
            if re.match(r'(\s*\w+\s*)=(\s*\(?[\-\+]?\w+[\+\-\*\/]?)',line):
                #mathematical exression sentence, such as: L1=2.1; angle=pi/24
                line_type = 'MATH_EXPR'
            elif re.match(r'^\s*%',line):
                #rpn sto expression
                line_type = 'MATH_EXPR_RPN_STO'
            elif re.match(r'(\s*\w+\s*):(\s*LINE\s*)=',line,re.IGNORECASE):
                # print(line)
                line_type = 'LINE'
            elif re.match(r'(\s*\w+\s*):(\s*[a-zA-Z]+\s*)', line):
                # print(line)
                line_type = 'ELEMENT'
            else:
                print('Unrecognized MATH_EXPR, ELEMENT, LINE name: ',line)
                sys.exit()

            if line_type == 'MATH_EXPR':
                # use exec to run the expression
                line = line.replace(' ','')
                exec(line.lower())   #use lower() in case pi is used.
            
            elif line_type == 'MATH_EXPR_RPN_STO':
                #like: % 1 2 / sto Lq
                line = line.lower()            # variables use lower case
                line = line.replace('%','')       
                
                # rpn calculation section
                # ---------------------------------
                # replace: mult, ...
                # OTHER CASES NOT INCLUDED YET
                line = line.replace(' mult ',' * ')                    
                
                stack = []
                for c in line.split():
                    if c in '+-*/':
                        i2 = stack.pop()
                        i1 = stack.pop()
                        #print(i1,c,i2)
                        #print( eval('%s'*3 % (i1,c,i2)) )
                        stack.append(eval('%s'*3 % (i1,c,i2)))
                    else:
                        stack.append(c)
                # ---------------------------------
                
                #finally, Lq=0.5                   
                exec("%s=%s" %(stack[-1],stack[0]))
                
            elif line_type == 'ELEMENT': 
                # handle with element line    
                tmp = line.split(',')
                
                # get element name and type
                elem_name = tmp[0].split(':')[0].upper()  #upper case for element name
                elem_type = tmp[0].split(':')[1].upper()  #upper case all elem_type
                
                # delete white space
                elem_name = elem_name.replace(' ','')
                elem_type = elem_type.replace(' ','')
                
                elem = dict()           
                elem['NAME'] = elem_name
                elem['TYPE'] = elem_type
                
                for para in tmp[1:]:
                    para_name = para.split('=')[0].upper().replace(' ','')  #use upper case for all para_name                                        
                    para_value = para.split('=')[1].lower()   #variables use lower case
                                        
                    # remove " or ' 
                    para_value = para_value.replace('"','')   
                    para_value = para_value.replace('\'','')
                    
                    # be careful of such as: coord_conv="NORMAL", if happened
                    # normal=0.2 in math expression, then 0.2 will be replaced.
                    try:
                        # normal math expression case: 
                        eval(para_value.lower())  #lower, in case pi is used.                  
                    except:
                        # rpn expression case
                        # ---------------------------------
                        # replace: mult, ...
                        # OTHER CASES NOT INCLUDED YET
                        para_value = para_value.replace(' mult ',' * ')                    
                        
                        stack = []
                        for c in para_value.split():
                            if c in '+-*/':
                                i2 = stack.pop()
                                i1 = stack.pop()
                                #print(i1,c,i2)
                                #print( eval('%s'*3 % (i1,c,i2)) )
                                stack.append(eval('%s'*3 % (i1,c,i2)))
                            else:
                                stack.append(c)
                        # ---------------------------------
                        para_value = stack[0]
                        
                    else:
                        para_value = eval(para_value.lower())
                                   
                    # back to str type
                    # for no-math expression, such as: option="zdE", back to upper
                    elem[para_name] = str(para_value).upper().replace(' ','')
                    
                lattice[elem_name] = elem   
                
            elif line_type == 'LINE':
                beamline = dict()
                line = line.replace(' ','')
                
                beamline_name = line.split(':')[0].upper() #upper case for beamline name
                beamline_line = line.split(':')[1].split('=')[1]
                beamline_line = beamline_line.strip('(').strip(')') #delete ()
                
                beamline['NAME'] = beamline_name
                beamline['TYPE'] = 'LINE'
                beamline['LINE'] = beamline_line.upper()
                
                lattice[beamline_name] = beamline
                        
        return(lattice) 

   
    def _expandline(self, line, lines):
        '''
        expand nested N*FODO to (FODO, FODO,...)
        expand FODO to (QF1,D1,QD1,....)
        '''  
        # lines = self.get_brieflines()
        
        lattice = self._get_lattice(lines)
        
        out = []        
        for item in line:
            if re.match(r'^\d+\*',item):  #4*FODO case
                tmp = item.split('*')
                tmp = int(tmp[0]) * [tmp[1]]                
                out.extend(self._expandline(tmp,lines))
            
            elif re.match(r'^[a-zA-Z]+\*',item):   #FODO*4 case
                tmp = item.split('*')
                tmp = int( tmp[1]) * [tmp[0]]
                out.extend( self._expandline(tmp,lines))
                print('ATTENTION: If in Elegant, you should use N*FODO, not FODO*N in .lte file.')    
            
            elif lattice[item]['TYPE'] == 'LINE':
                item = lattice[item]['LINE'].split(',')
                out.extend( self._expandline(item,lines)) # in case item is still LINE, recursion
                
            else:
                out.append(item)                    
        return out   
    
    def _is_number(self, s:str):
        '''
        To judge whether input string s is a number or not.
        '''
        try:
            float(s)
            return True
        except ValueError:
            return False      
        
    def rpn_cal(self, rpn_str:str):
        '''
        calculate the rpn expression:  
        case 1:
        rpn_str='1 3 /'
        return is [0.333]
        case 2:
        rpn_str="% 1 3 / sto abc"
        return is [0.333, 'sto', 'abc']        
        '''
        # -------------------------------------------
        # replace 'mult' with *
        # OTHER rpn expressions are not included YET
        # -------------------------------------------
        rpn_str = rpn_str.replace(' mult ',' * ')        
        stack = []
        for c in rpn_str.split():
            if c in '+-*/':
                i2 = stack.pop()
                i1 = stack.pop()
                #print(i1,c,i2)
                #print( eval('%s'*3 % (i1,c,i2)) )
                stack.append(eval('%s'*3 % (i1,c,i2)))
            else:
                stack.append(c)
        return stack

if __name__=='__main__':

    # debug 
    # ========
    file_name = './debugExample/lattice.lte'
    line_name = 'line'

    lte = lattice_parser(file_name,line_name)
    lines = lte.get_brieflines()
    trackline = lte.get_trackline(lines)
