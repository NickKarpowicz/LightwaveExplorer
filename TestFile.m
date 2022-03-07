fid = fopen('TestFile_ExtIn.dat','rb'); 
TestFile_ExtIn = fread(fid,  65536, 'double'); 
TestFile_ExtIn = reshape(TestFile_ExtIn,[256 128 2]); 
fclose(fid); 
fid = fopen('TestFile_ExtOut.dat','rb'); 
TestFile_ExtOut = fread(fid,  65536, 'double'); 
TestFile_ExtOut = reshape(TestFile_ExtOut,[256 128 2]); 
fclose(fid); 
