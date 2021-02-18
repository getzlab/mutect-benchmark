# Google VM + VNC + IGV

## Google VM

The pub key email username should match the username in the server!



## Run VNC

On Google VM,

```bash
vncserver :7 -geometry 1200x800 -depth 24 -rfbauth ~/.vnc/passd
```

The number after the colon, plus `5900` (in this case, 5907) is the port we specified!



On local computer,

```
ssh -L 5907:localhost:5907 -N -f -l qing EXTERNAL_IP
```

then run `vncviewer`and type the passwd to proceed on Ubuntu, `open vnc://localhost:5907` on Mac.



**Sidenote** When the passwd file is missing for some reason - repopulate it.

```bash
echo MYVNCPASSWORD | vncpasswd -f > ~/.vnc/passd
```



## IGV

note that IGV needs java 11. (for instructions [here](<https://www.linuxuprising.com/2019/06/new-oracle-java-11-installer-for-ubuntu.html>), JDK must be downloaded from Oracle website)

In the VNC terminal - 

```
java --module-path=lib -Xmx4G @igv.args --module=org.igv/org.broad.igv.ui.Main
```

create an alias will help - remember first to locate into the IGV_Linux directory



# misc

* Enable copy and paste in VNC ([link](<https://superuser.com/questions/1081489/how-to-enable-text-copy-and-paste-for-vnc>))
* Use google bucket



