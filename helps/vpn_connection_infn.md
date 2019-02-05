# Import a VPN connection in Linux for INFN connection
These instructions refer mainly to Naples Infn connection, which uses two files, one configuration file and one certificate file.
To import a VPN connection, we need the OpenVpn plugin to the Network Manager. The following step should apply:
* sudo apt-get install openvpn
* sudo apt-get install network-manager-openvpn-gnome
* It might be needed to restart the network manager with sudo /etc/init.d/networking restart
* Launch the Network manager and go to VPN connections
* Import VPN from file
* Select the .ovpn file
* Insert the certificate file when needed