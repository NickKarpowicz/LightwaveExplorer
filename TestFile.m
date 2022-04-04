fid = fopen('TestFile_ExtIn.dat','rb'); 
TestFile_ExtIn = fread(fid, 327680, 'double'); 
TestFile_ExtIn = reshape(TestFile_ExtIn,[512 320 2]); 
fclose(fid); 
fid = fopen('TestFile_ExtOut.dat','rb'); 
TestFile_ExtOut = fread(fid, 327680, 'double'); 
TestFile_ExtOut = reshape(TestFile_ExtOut,[512 320 2]); 
fclose(fid); 
