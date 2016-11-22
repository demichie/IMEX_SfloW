clear all;

[filename, pathname] = uigetfile( '*.p*','Pick a file','MultiSelect','on');

if ( ~iscell(filename) )
    
    nfiles = 1;
    
else
    
    nfiles = length(filename);

end

figure

for i=1:nfiles,
    
    if ( nfiles == 1 )
        
        if ( ismac)
            
            importfile_mac(strcat(pathname,filename));
            
        else
            
            importfile_colima(strcat(pathname,filename));
            
        end
        
    else
        
        if ( ismac)
            
            importfile_mac(strcat(pathname,filename{i}));
            
        else
            
            importfile_colima(strcat(pathname,filename{i}));
            
        end
        
    end
    
    a = (reshape(data',5,[]))';
    
    p = a(:,1:5);
    
    x0 = str2double(strrep(textdata{1},'x0',''));
    dx = str2double(strrep(textdata{2},'dx',''));
    cells = str2double(strrep(textdata{3},'cells',''));
    t = str2double(strrep(textdata{4},'t',''));
    
    N = size(p,1);
    
    x = x0+dx* (0.5 + 0:N);
        
    
    %% Plot density of the first phase
    
    subplot(nfiles,2,2*(i-1)+1);
    
    plot(x,p(1:N,4),x,p(1:N,5));
    
    hold all;
    ylabel('Height (m)');
    xlabel('Position (m)');
   
    str_title = ' Height --- t=';
    title([str_title,num2str(t)]);
    
    subplot(nfiles,2,2*(i-1)+2);
    
    plot(x,p(1:N,3));
    str_title = ' Velocity --- t=';
    title([str_title,num2str(t)]);
    ylabel('Velocity (m)');
    xlabel('Position (m)');
    
end


