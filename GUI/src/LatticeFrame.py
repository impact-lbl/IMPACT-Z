'''
'''
import tkinter as tk
import sys, threading
_TextWidth = 70
_TextHeight = 9
ELEMENT_TYPE    = {'DRIFT'      :   '0',
                   'QUADM'      :   '1',
                   'CONSTFOC'   :   '2',
                   'SOLENOID'   :   '3',
                   'DIPOLE'     :   '4',
                   'MULTIPOLE'  :   '5',
                   'DTL'        :   '101',
                   'CCDTL'      :   '102',
                   'CCL'        :   '103',
                   'SC'         :   '104',
                   'SOLRF'      :   '105',
                   'EMFLD'      :   '110',
                   'EMFLDCART'  :   '111',
                   'EMFLDCYL'   :   '112',
                   'EMFLDANA'   :   '113' }

ELEMENT_COLOR    = {'DRIFT'     :   '#444444',      #SILVER
                   'QUADM'      :   '#FFDC00',      #YELLOW
                   'CONSTFOC'   :   '#FF851B',      #ORANGE
                   'SOLENOID'   :   '#2ECC40',      #GREEN
                   'DIPOLE'     :   '#0074D9',      #BLUE
                   'MULTIPOLE'  :   '#001f3f',      #NAVY
                   'DTL'        :   '#AAAAAA',      #GRAY
                   'CCDTL'      :   '#39CCCC',      #TEAL
                   'CCL'        :   '#3D9970',      #OLIVE
                   'SC'         :   '#7FDBFF',      #AQUA
                   'SOLRF'      :   '#01FF70',      #LIME
                   'EMFLD'      :   '#FF4136',      #RED
                   'EMFLDCART'  :   '#B10DC9',      #PURPLE
                   'EMFLDCYL'   :   '#F012BE',      #FUCHSIA
                   'EMFLDANA'   :   '#85144b' }     #MAROON


ELEMENT_TYPE_INVERSE = dict(map(lambda t:(t[1],t[0]), ELEMENT_TYPE.items()))

class LatticeFrame(tk.Frame):
    """Lattice"""

    def __init__(self, master=None, cnf={}, **kw):
        '''See the __init__ for tk.Frame for most of this stuff.'''
        tk.Frame.__init__(self, master, cnf, **kw)
        self.createWidgets()  

    def createWidgets(self): 
        self.lattice_sv = tk.Scrollbar(self, orient=tk.VERTICAL)
        self.lattice_sh = tk.Scrollbar(self, orient=tk.HORIZONTAL)
  
        self.latticeText = tk.Text(self,
                               width = _TextWidth,height=_TextHeight,
                               wrap='none',
                               yscrollcommand=self.lattice_sv.set,
                               xscrollcommand=self.lattice_sh.set)

        self.lattice_sv.config(command=self.latticeText.yview)  
        self.lattice_sh.config(command=self.latticeText.xview)  
  
        self.lattice_sv.pack(fill="y", expand=0, side=tk.RIGHT, anchor=tk.N)  
        self.lattice_sh.pack(fill="x", expand=1, side=tk.BOTTOM, anchor=tk.N)  
        self.latticeText.pack(fill="both", expand=1, side=tk.TOP)  
  

        self.latticeText.bind("<Control-Key-a>", self.selectText)  
        self.latticeText.bind("<Control-Key-A>", self.selectText)  
        self.latticeText.insert(tk.END,'20.0d0 1 1 0 0.d0 0.12 0.12 /')

    
    def get(self, index1, index2=None):
        return self.latticeText.get(index1, index2)
        
    def selectText(self, event):
        self.latticeText.tag_add(tk.SEL, "1.0", tk.END)  
        return 'break'

class LatticeFrameC(tk.Frame):
    """Lattice"""
    class TextLineNumbers(tk.Canvas):
        def __init__(self, *args, **kwargs):
            tk.Canvas.__init__(self, *args, **kwargs)
            self.textwidget = None
    
        def attach(self, text_widget):
            self.textwidget = text_widget
    
        def redraw(self, *args):
            '''redraw line numbers'''
            self.delete("all")
    
            i = self.textwidget.index("@0,0")
            while True:
                dline= self.textwidget.dlineinfo(i)
                if dline is None: break
                y = dline[1]+15
                linenum = str(i).split(".")[0]
                self.create_text(30,y,anchor="ne", text=linenum)
                i = self.textwidget.index("%s+1line" % i)
    
    class CustomText(tk.Text):
        def __init__(self, *args, **kwargs):
            tk.Text.__init__(self, *args, **kwargs)
    
            self.tk.eval('''
                proc widget_proxy {widget widget_command args} {
    
                    # call the real tk widget command with the real args
                    set result [uplevel [linsert $args 0 $widget_command]]
    
                    # generate the event for certain types of commands
                    if {([lindex $args 0] in {insert replace delete}) ||
                        ([lrange $args 0 2] == {mark set insert}) || 
                        ([lrange $args 0 1] == {xview moveto}) ||
                        ([lrange $args 0 1] == {xview scroll}) ||
                        ([lrange $args 0 1] == {yview moveto}) ||
                        ([lrange $args 0 1] == {yview scroll})} {
    
                        event generate  $widget <<Change>> -when tail
                    }
    
                    # return the result from the real widget command
                    return $result
                }
                ''')
            self.tk.eval('''
                rename {widget} _{widget}
                interp alias {{}} ::{widget} {{}} widget_proxy {widget} _{widget}
            '''.format(widget=str(self)))
            
    def __init__(self, master=None, cnf={}, **kw):
        '''See the __init__ for tk.Frame for most of this stuff.'''
        tk.Frame.__init__(self, master, cnf, **kw)
        self.createWidgets()  

    def createWidgets(self):   
        self.lattice_sv = tk.Scrollbar(self, orient=tk.VERTICAL)
        self.lattice_sh = tk.Scrollbar(self, orient=tk.HORIZONTAL)
  
        self.latticeText = self.CustomText(self,
                               width = _TextWidth,height=_TextHeight,
                               wrap='none',
                               yscrollcommand=self.lattice_sv.set,
                               xscrollcommand=self.lattice_sh.set)
        
        self.latticeText.tag_configure("bigfont", font=("Helvetica", "24", "bold"))
        self.linenumbers = self.TextLineNumbers(self, width=30,height=_TextHeight,)
        self.linenumbers.attach(self.latticeText)
        
        self.lattice_sv.config(command=self.latticeText.yview)  
        self.lattice_sh.config(command=self.latticeText.xview)
        
        self.title = tk.Text(self,width = _TextWidth,height=1,
                               wrap='none',bg='#ededed',borderwidth=0)
        self.titleT()

        self.lattice_sv.pack(fill="y", expand=0, side=tk.RIGHT, anchor=tk.N)
        self.linenumbers.pack(fill="y",expand=1, side=tk.LEFT,  anchor=tk.N)
        self.lattice_sh.pack(fill="x", expand=1, side=tk.BOTTOM,anchor=tk.N)
        self.title.pack(fill="both", expand=1, side=tk.TOP) 
        self.latticeText.pack(fill="both", expand=1, side=tk.TOP)  
  

        self.latticeText.bind("<Control-Key-a>", self.selectText)  
        self.latticeText.bind("<Control-Key-A>", self.selectText)
        self.latticeText.bind("<<Change>>", self._on_change)
        self.latticeText.bind("<Configure>", self._on_change)
    
        self.latticeTextHide = tk.Text(self,
                               width = _TextWidth,height=3,
                               wrap='none')
        #self.latticeTextHide.pack(fill="both", expand=1, side=tk.TOP)

        self.latticeTextHide.insert(tk.END,'20.0 1 1 0 0.d0 0.12 0.12 /')
        self.update()
        
        for ele in ELEMENT_COLOR.keys():
            self.latticeText.tag_config(ele, foreground=ELEMENT_COLOR[ele])
        self._on_change('change')
        
    def _on_change(self, event):
        self.linenumbers.redraw()
        for ele in ELEMENT_COLOR.keys():
            self.search(self.latticeText,  ele,   ele)
    def get(self, index1, index2=None):
        return self.latticeText.get(index1, index2)
    
    def getHide(self):
        self.updateHide()
        return self.latticeTextHide.get('0.0', tk.END)
    
    def update(self):
        self.latticeText.delete('0.0', tk.END)
        text = self.latticeTextHide.get('1.0', tk.END).splitlines()
        for line in text:
            if line.strip()!='':
                a=self.convertNtoW(line)
                self.latticeText.insert('end',a)
            
    def updateHide(self):
        self.latticeTextHide.delete('0.0', tk.END)
        text = self.latticeText.get('1.0', tk.END).splitlines()
        for line in text:
            if line.strip()!='':
                a=self.convertWtoN(line)
                self.latticeTextHide.insert('end',a)
       
    def selectText(self, event):
        self.latticeText.tag_add(tk.SEL, "1.0", tk.END)  
        return 'break'
    
    def convertNtoW(self,s1):
        strSet = s1.split()
        try:
            eleName = ELEMENT_TYPE_INVERSE[strSet[3]]
            wordFormat = '{:<10}'.format(eleName)
            wordFormat +=  ' ' + '{:<10}'.format(str(float(strSet[0])))
            for s in strSet[1:3]:
                wordFormat +=  ' ' + '{:<10}'.format(s)
            for s in strSet[4:]:
                wordFormat +=  ' ' + s
            return wordFormat+'\n'
        except:
            return s1+'\n'
            
    
    def convertWtoN(self,s1):
        strSet = s1.split()
        try:
            numberFormat = strSet[1] + ' ' + strSet[2] + ' ' + strSet[3] + ' ' + str(ELEMENT_TYPE[strSet[0]])
            for s in strSet[4:]:
                numberFormat +=  ' ' + s
            return numberFormat+'\n'
        except:
            return s1+'\n'
        
    def titleT(self):
        self.title.config(state='normal')
        wordFormat = '{:<10}'.format('Name')
        wordFormat +=  ' ' + '{:<10}'.format('Length')
        wordFormat +=  ' ' + '{:<10}'.format('')
        wordFormat +=  ' ' + '{:<10}'.format('')
        wordFormat +=  ' ' + '{:<15}'.format('V1,V2,V3...')
        self.title.delete('0.0', tk.END)
        self.title.insert('0.0',wordFormat)
        self.title.config(state='disabled')
    def titleZ(self):
        self.title.config(state='normal')
        wordFormat = '{:<10}'.format('Name')
        wordFormat +=  ' ' + '{:<10}'.format('Length')
        wordFormat +=  ' ' + '{:<10}'.format('Step')
        wordFormat +=  ' ' + '{:<10}'.format('Map Step')
        wordFormat +=  ' ' + '{:<15}'.format('V1,V2,V3...')
        self.title.delete('0.0', tk.END)
        self.title.insert('0.0',wordFormat)
        self.title.config(state='disabled')
    
    def search(self, text_widget, keyword, tag):
        pos = '1.0'
        while True:
            idx = text_widget.search(keyword, pos, tk.END)
            if not idx:
                break
            pos = '{}+{}c'.format(idx, len(keyword))
            text_widget.tag_add(tag, idx, pos)
            
class ConsoleText(tk.Text):
    '''A Tkinter Text widget that provides a scrolling display of console stdout.'''

    class IORedirector(object):
        '''A general class for redirecting I/O to this Text widget.'''
        def __init__(self,text_area):
            self.text_area = text_area

    class StdoutRedirector(IORedirector):
        '''A class for redirecting stdout to this Text widget.'''
        def write(self,mystr):
            self.text_area.write(mystr,False)

    class StderrRedirector(IORedirector):
        '''A class for redirecting stderr to this Text widget.'''
        def write(self,mystr):
            self.text_area.write(mystr,True)

    def __init__(self, master=None, cnf={}, **kw):
        '''See the __init__ for Tkinter.Text for most of this stuff.'''

        tk.Text.__init__(self, master, cnf, **kw)

        self.started = False
        self.write_lock = threading.RLock()

        #self.tag_configure('STDOUT',background='black',foreground='white')
        #self.tag_configure('STDERR',background='black',foreground='red')
        
        #self.console_sv = tk.Scrollbar(master, orient=tk.VERTICAL)
        #self.console_sv.config(command=self.yview)
        #self.console_sv.pack(fill="y", expand=0, side=tk.RIGHT, anchor=tk.N) 
        self.width = _TextWidth
        self.height=_TextHeight,
        self.config(state=tk.NORMAL)
        self.bind("<1>", lambda event: self.focus_set())
        #self.bind('<Key>',lambda e: 'break') #ignore all key presses

         
    def start(self):

        if self.started:
            return

        self.started = True

        self.original_stdout = sys.stdout
        #self.original_stderr = sys.stderr

        stdout_redirector = ConsoleText.StdoutRedirector(self)
        #stderr_redirector = ConsoleText.StderrRedirector(self)

        sys.stdout = stdout_redirector
        #sys.stderr = stderr_redirector

    def stop(self):

        if not self.started:
            return

        self.started = False

        sys.stdout = self.original_stdout
        #sys.stderr = self.original_stderr

    def write(self,val,is_stderr=False):

        #Fun Fact:  The way Tkinter Text objects work is that if they're disabled,
        #you can't write into them AT ALL (via the GUI or programatically).  Since we want them
        #disabled for the user, we have to set them to NORMAL (a.k.a. ENABLED), write to them,
        #then set their state back to DISABLED.

        self.write_lock.acquire()

        self.update_idletasks()
        self.insert('end',val)
        self.see('end')

        self.write_lock.release()
    
