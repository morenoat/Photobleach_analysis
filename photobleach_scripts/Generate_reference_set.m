function Generate_reference_set(output_name_ID_molecules,ID_molecules,plot_each_image,last_round_plot,path_to_data)


point_distance=[4;3;2;1;0.5;0.5;0.5;0.5];
start_input=[output_name_ID_molecules '_ID'];
name_ref_out=output_name_ID_molecules;
d=uigetdir('','Select Input-folder'); %select the input-folder that contains the subfolders
cd(d);
folder_path=pwd;
rename_image_files
file_list=dir(fullfile(folder_path, '*.tif'));
number_of_images=numel((file_list));
if ID_molecules==true
    Identify_reference_molecules(number_of_images,20,start_input,10,6);
end
cd([pwd '/' start_input]);
for jk=1:size(point_distance,1)
    if jk==1
        refinp=[start_input '.mat'];
        refout=name_ref_out;
    else
        refinp=[name_ref_out '.mat'];
        refout=name_ref_out;
    end
    if plot_each_image==false && jk==size(point_distance,1)
        calibration_refine(refinp,number_of_images,8,15,point_distance(jk),plot_each_image,refout,path_to_data,'Last_cycle',true,'plot_last_batch',last_round_plot);
    elseif plot_each_image==true
        calibration_refine(refinp,number_of_images,8,15,point_distance(jk),plot_each_image,refout,path_to_data,'Last_cycle',true);
    else
        calibration_refine(refinp,number_of_images,8,15,point_distance(jk),plot_each_image,refout,path_to_data);
    end

end
