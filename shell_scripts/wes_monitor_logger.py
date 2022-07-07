#Written by jacobg 2022
import requests
import sys
import os

#GET (my) instance_name and the wes_monitor server IP and port from bash envars
INSTANCE_NAME=os.environ.get("HOSTNAME")
SERVER_IP_PORT=os.environ.get("WES_MONITOR_IP_PORT")

def log_handler(msg):
    level = msg["level"]
        
    #RUNS AT RUN START
    if level == "run_info":
        #msg['msg'] is a table, whose last line is either
        #OLDER version of snakemake: last line = total number of jobs
        #NEWER version- a row of the table, e.g. total   6     1     1
        # where the 1st col is the total
        last_line = msg['msg'].split("\n")[-1]
        if len(last_line.split()) == 1:
            #OLDER snakemake
            total = last_line.strip()
        else:
            #NEWER snakemake
            #total= msg['msg'].strip().split("\n")[-1].split()[1]
            total= last_line.split()[1]
            
        #NOTE: need to -1 from the total to correct for the "all" job
        total = str(int(total) - 1)
        #print("we have %s total jobs" % total)
        try:
            r=requests.put("http://%s/update/%s" % (SERVER_IP_PORT,INSTANCE_NAME), data={'num_steps':total})
            #print(r.content)
        except:
            print("WES Monitor server down %s" % SERVER_IP_PORT)

        
    #RUNS WHEN JOBS COMPLETE, INCLUDING RUN COMPLETION 
    elif level == "progress":
        step_count = msg["done"]
        if step_count == msg["total"]:
            #print("we did it!")
            try:
                r=requests.put("http://%s/update/%s" % (SERVER_IP_PORT,INSTANCE_NAME), data={'status':"COMPLETE"})
            except:
                print("WES Monitor server down %s" % SERVER_IP_PORT)
        else:
            #print("still going")
            #print(msg["done"])
            try:
                r=requests.put("http://%s/update/%s" % (SERVER_IP_PORT,INSTANCE_NAME), data={'step_count':step_count})
            except:
                print("WES Monitor server down %s" % SERVER_IP_PORT) 

    elif level == "error":
        #extract which rule has caused the error
        #relay info to wes monitor
        #print("uh oh! we have an error!")
        try:
            r=requests.put("http://%s/update/%s" % (SERVER_IP_PORT,INSTANCE_NAME), data={'status':"ERROR"})
        except:
            print("WES Monitor server down %s" % SERVER_IP_PORT) 
