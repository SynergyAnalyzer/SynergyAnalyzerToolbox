# Data Requirements

## Data Structure

The toolbox assumes data are organized in a specific structure. This is an example of the structure for a reaching paradigm.
Data for the demo are from the study "Mixed matrix factorization: a novel algorithm for the extraction of kinematic-muscular synergies" by Alessandro Scano and collaborators, published in 2022 on Journal of Neurophysiology.
Here we selected a data subset of one participant, reaching towards 8 targets on the sagittal plane. There are 10 repetitions of each trial.

Loading raw EMG and kinematic data (reaching to 8 target in sagittal plane)

```matlab
load('data.mat');
disp(data)
```

  1×80 struct array with fields:

    emg
    pos
    info
    emgtime
    postime

The ```emg``` field has ```n * t``` dimensions where ```n``` is the number of channels and ```t``` the number of time samples. Therefore the ```emgtime``` field has ```t``` samples. 
Similarly the ```postime``` field has ```k``` samples and pos has ```m * k``` elements, where ```m``` is the number of channels.
The ```info``` field contains the information about the protocol and how the data have been collected. For example:

```matlab
disp(data(1).info)
```

      plane: 1
      start: 0
     target: 1
    t_onset: 0
      t_end: 0.5100

```plane``` indicates the plane in which the targets were displayed (1->frontal, 2-> sagittal).
```start``` indicates the target from which the movement started.
```target``` indicates the target to reach.
```t_onset``` indicates the time of the onset of the movement, (i.e. 10% of max speed).
```t_end``` indicates the time of the offset of the movement, (i.e. 10% of max speed after the peak).

NB: the software assumes that these fields are completed in the data structure.


## Time alignment                

In addition EMG and kinematic data need to be aligned in time, preferably to the onset time, i.e. t_onset = 0.


## Trial duration

In order to subtract the tonic component to the EMG signal and extract the phasic activity, trial duration should include 400ms before onset time and 400ms after the end of the movement.

# From data struct to info struct

SynergyAnalyzer class has a property data and a property info and both these properties are transformed from the structures in input by the SynergyAnalyzer class constructor. To this end we include a function that help in this transformation.
The ```getInfo``` function transform the info field of the data structure into the info struct as required by the SynergyAnalyzer property.
The output of the ```getInfo``` function needs to have the following fields for each trial ```i``` in the data structure:

```matlab
    info(i).id          %this is the trial number
    info(i).type        %this could be a vector of n different types
    info(i).selected    %this indicates if the trial could be included in the analysis
   
    info(i).events.code = [13 14];                                      % these are the codes for the time events (we included two but more can be added)
    info(i).events.time = [data(i).info.t_onset data(i).info.t_end];    %these are the times at which the events occurred
    
```

In our example the ```getInfo``` function set as follows:

```matlab

ntrial = length(data);

for i=1:ntrial
    % id
    info(i).id = i;
    
    % type
    info(i).type(1) = data(i).info.plane;   % plane type
    
    if data(i).info.start==0   % center-out/out-center movement, direction of motion
        info(i).type(2) = 1;   
        info(i).type(3) = data(i).info.target;
    else
        info(i).type(2) = 0;
        info(i).type(3) = mod(data(i).info.start+3,8)+1;  % assumes 8 targets
    end
    
    % selected
    info(i).selected = 1;
    
    % events
    info(i).events.code = [13 14];
    info(i).events.time = [data(i).info.t_onset data(i).info.t_end];
    
    
end

```