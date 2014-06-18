#temporal difference RT model
#try on one run
library(fitclock)

exSubj <- clockdata_subject(subject_ID="008", dataset=clocksubject_fMRI_008jh)
run1 <- exSubj$runs[[1L]]

rts <- run1$RTraw
rewards <- run1$Reward

#adapted ludvig code
#Simple Code for CSC/PR/MS TD Model with Actions.

#Model Parameters
numms = 12         #Number of microstimuli per stimulus
numstimuli = 1     #Number of stimuli, including the reward/US (here we just have US)
alpha = 0.05       #Step-size
#decay = 0.985      #Memory Trace decay rate
decay = 0.999      #Memory Trace decay rate
gamma = 0.97       #Discount factor
lambda = 0.95      #Eligibility trace decay rate
sigma = 0.08       #Microstimulus width
theta = 0.25       #Action Threshold
upsilon = 0.9      #Action Decay
presence = 0.2     #Presence microstimulus level

# Experiment Code:
# CS2 is used for Overshadowing and Blocking

numtrials = length(rts)    
numtimesteps = 4000 #1ms bins in a 4000ms trial

## Initialize Data Vectors:
x = rep(0, (numms + 1)*numstimuli)         #Stimulus representation
w = rep(0, (numms + 1)*numstimuli)         #Weight vector
e = rep(0, (numms + 1)*numstimuli)         #Eligibility Traces
delta = matrix(data=0, nrow=numtrials, ncol=numtimesteps)    #TD Errors
value = matrix(data=0, nrow=numtrials, ncol=numtimesteps)    #Value functions
action = matrix(data=0, nrow=numtrials, ncol=numtimesteps)   #Response Levels
ms = matrix(data=0, nrow=numtimesteps, ncol=numms)           #Microstimulus levels
#maxaction = zeros (1, numtrials);           

stimulustime <- matrix(rts, ncol=1) #stimulus onset time, here just the US

## Microstimulus Generation (NIPS2009 Version)
trace = 1;
for (timestep in 1:numtimesteps) {
  for (micro in 1:numms) {
    ms[timestep, micro] = trace * ((1/sqrt(2*pi))*exp(-((trace-(micro/numms))^2)/(2*(sigma^2))));    
  }
  trace = trace * decay;
}

#plot microstimuli
library(reshape2)
library(ggplot2)

ms_melt <- melt(ms, varnames=c("timestep", "msnumber"))
ms_melt$msnumber <- factor(ms_melt$msnumber)
ggplot(ms_melt, aes(x=timestep, y=value, color=msnumber)) + geom_line()

##fit data
##Run Simuluation:
  
for (trial in 1:numtrials) {
  oldvalue = 0;               # Reset on every trial
  #oldreward = 0;
  oldaction = 0;
  e = rep(0, (numms + 1)*numstimuli)
  rewardtime = rts[trial];
  
  for (timestep in 1:numtimesteps) {
    if (timestep == rewardtime) {
      reward = rewards[trial]
    } else {    
      reward = 0
    }
    
    #microstimulus representation of stimulus
    for (stimulus in 1:numstimuli) {
      
      if (timestep >= stimulustime[trial, stimulus]) {
        #noise = rand*0.02 -.01; (unused)
        noise = 0;
        x [(stimulus-1)*numms +1:(stimulus-1)*numms +numms] = ms[timestep-stimulustime[trial, stimulus]+1,] + noise;     
      } else {
        x [(stimulus-1)*numms +1:(stimulus-1)*numms +numms] = 0;
        x [numstimuli * numms + stimulus] = 0;
      }
    }

    value (trial, timestep) = crossprod(x,w)
    
    #Action Selection:
    action[trial, timestep] = upsilon * oldaction + max (0, oldvalue - theta)
    oldaction=action[trial, timestep]
    
    #Learning Algorithm:
    delta[trial, timestep] = reward + (gamma * value (trial, timestep)) - oldvalue #TD Learning
    
    w = w + (alpha * delta[trial, timestep] * e)
    e = x + (gamma * lambda * e)
    oldvalue = crossprod(x,w) #Or oldvalue = value. Difference in which weights are used.
    #oldreward = reward
  
  }

  
      
  
}
    

  

#%% Real-time plotting:
#    %     maxaction(k) = max(max(action));
#    plot(value (trial, :), 'r');
#hold on;
#plot(delta (trial, :), 'b');
#plot(action (trial, :), 'g');
#hold off;
#drawnow;
#end
