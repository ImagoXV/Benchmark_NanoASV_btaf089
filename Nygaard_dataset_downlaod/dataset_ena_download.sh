#/bin/bash

nygaard=(SAMEA4883696 SAMEA4883697 SAMEA4883698 SAMEA4883699 SAMEA4883700 SAMEA4883701 SAMEA4883702 SAMEA4883703 SAMEA4883704 SAMEA4883705 SAMEA4883706)

#For this script to work, you need to adapt the path to the enBrowserTools binary

for ena in "${nygaard[@]}"; do
	/path/to/enaBrowserTools-1.7.1/python3/enaDataGet $ena
done

