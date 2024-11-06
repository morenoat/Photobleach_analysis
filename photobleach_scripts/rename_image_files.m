%%%This must be used in the directory with images

a =pwd;
A =dir( fullfile(a, '*.tif') );
fileNames = { A.name };
pattern='ng';
for iFile = 1 : numel( A )
    if isempty(regexp(fileNames{ iFile }, pattern, 'once'))
        newName = fullfile(a, sprintf( 'ng%d.tif', iFile ) );
        movefile( fullfile(a, fileNames{ iFile }), newName );
    end
end