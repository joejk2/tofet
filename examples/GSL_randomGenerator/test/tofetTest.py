import unittest, os
from subprocess import Popen, PIPE, call 
import re

class tofetTest(unittest.TestCase):
    
    def setUp(self):
        self.test_output = ""
    
    def tearDown(self):
        os.system("rm -f tofetTest.out")
        pass

    def testToF(self):
        """Check simple ToF simulations"""

        self.sim =  "../tof/tof.sim"
        self.ref =  "../tof/tof_0.out"
        self.edge = "../tof/scl.edge"
        self.xyz  = "../tof/scl.xyz"

        self.read_reference_file()
        self.run_simulations()
        assert self.compare_values("> MOBILITY FROM COLLECTION TIMES (cm^2/V.s)=")
        assert self.compare_values("> PROBABILITY OF HOPPER BEING COLLECTED DURING RUN =")

    def testRegenerateOcc(self):
        """Check regenerate mode, with occupation probabilities"""

        self.sim =  "../regenerate_occ/regenerate_occ.sim"
        self.ref =  "../regenerate_occ/no_traps.out"
        self.edge = "../regenerate_occ/scl_no_traps.edge"
        self.xyz  = "../regenerate_occ/scl.xyz"

        self.read_reference_file()
        self.run_simulations("tofetOccupation")

        assert self.compare_values("> MOBILITY FROM COLLECTION TIMES (cm^2/V.s)=")
        assert self.compare_values("> MOBILITY FROM TOTAL DISPLACEMENT AND TOTAL TIME (cm^2/V.s)=")
        assert self.compare_values("> PROBABILITY OF HOPPER BEING COLLECTED DURING RUN =")
        assert self.compare_series("> TOTAL OCCUPATION TIMES AND TIMES VISITED")
        #assert self.compare_files("../regenerate_occ/occProb.out", "occProb.out")

    def testFeT(self):
        """Check simple FeT simulations"""

        self.sim =  "../fet/fet.sim"
        self.ref =  "../fet/fet_0.out"
        self.edge = "../fet/scl_fet.edge"
        self.xyz  = "../fet/scl_fet.xyz"

        self.read_reference_file()
        self.run_simulations()

        assert self.compare_values("> TIME =")
        assert self.compare_values("> NUMBER OF HOPPERS LEFT =")
        assert self.compare_values( "> TOTAL NUMBER OF HOPPERS COLLECTED AT DRAIN =" )
        assert self.compare_values( "> TOTAL NUMBER OF HOPPERS INJECTED AT SOURCE =")
        assert self.compare_series("> OCCUPIED MOLECULES AT END OF SIMULATION")

    def run_simulations(self, binary = "tofet"):
        out = Popen( (binary, self.edge, self.xyz, self.sim), stdout=PIPE).stdout
        self.test_output = [line.strip() for line in out.readlines()]
        out.close()

    @staticmethod
    def get_value(outputSource_, stringList_, lineBegin_):
        for line in stringList_:
        #Find the line that contains the (escaped) string lineBegin_... :
            if re.match( re.escape(lineBegin_), line) is not None: #Needs to be escaped here...
                #..., strip lineBegin_ from the front, split at every " " 
                #... and return first remaining value 
                retval = line.strip(lineBegin_).split()[0]
                print outputSource_ + ": " + lineBegin_, retval
                return retval
        print outputSource_ + ":" + lineBegin_ + "Not Found!"
        return False

    def compare_series(self, lineBegin_):
        ref_index_begin, test_index_begin = "Not found" , "Not found"
        try:  
            ref_index_begin = self.ref_output.index(lineBegin_)   #Doesn't like being escaped here....
            test_index_begin = self.test_output.index(lineBegin_)
        except ValueError: 
            print "*** ERROR ***: Couldn't find", lineBegin_
            print "               ref_index, test_index: ", ref_index_begin + ",", test_index_begin
            return False
        
        list = []
        #Put everything in one list (kind of dirty)....
        for series in self.ref_output[ref_index_begin+1:], self.test_output[test_index_begin+1:]:
            for line in series:
                if ">" in line: break
                else: list.append(line)
        length = len(list)
        #Divide list into two...
        ref_series, test_series = [list[i:i + length/2] for i  in range(0, length, length/2)]
        
        return ref_series == test_series
 

    def compare_files(self, filename1_, filename2_):
        file1_, file2_ = open(filename1_), open(filename2_)
        data1_, data2_ = [line.strip() for line in file1_.readlines()], [line.strip() for line in file2_.readlines()] 
        file1_.close()
        file2_.close()
        
        if (data1_ == data2_):  return True
        else:                   return False

    def read_reference_file(self):
        ref_file = open(self.ref)
        self.ref_output = [line.strip() for line in ref_file.readlines()]
        ref_file.close()

    def compare_values(self, lineBegin_ ):

        refValue  = self.get_value(self.ref, self.ref_output, lineBegin_)
        testValue = self.get_value("test", self.test_output, lineBegin_)

        if (refValue == False) or (testValue == False):  
            print "refValue, testValue", refValue, ", ", testValue
            return False

        try:
            assert refValue == testValue 
        except AssertionError:
            print "*** ERROR *** : Values not identical!"
            print "                Line beginning: ", re.escape(lineBegin_)
            print "                Values: ", refValue + ", ", testValue  
            return False
        
        return True

if __name__ == '__main__':
    unittest.main()
