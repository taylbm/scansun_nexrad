import os
import sys
from dateutil import tz
from datetime import datetime, timedelta
import subprocess
import boto3
from suntime import Sun
from botocore import UNSIGNED
from botocore.config import Config

user_slice_start = int(sys.argv[1]) if len(sys.argv) > 1 else None
user_slice_end = int(sys.argv[2]) if len(sys.argv) > 2 else None


ROOT_PATH = '/home/btaylor/nexrad_data/'
NEXRAD_L2_BUCKET = 'noaa-nexrad-level2'
SITES =['KABR', 'KABX', 'KAKQ', 'KAMA', 'KAMX', 'KAPX',
        'KARX', 'KATX', 'KBBX', 'KBGM', 'KBHX', 'KBIS',
        'KBLX', 'KBMX', 'KBOX', 'KBRO', 'KBUF', 'KBYX',
        'KCAE', 'KCBW', 'KCBX', 'KCCX', 'KCLE', 'KCLX',
        'KCRP', 'KCXX', 'KCYS', 'KDAX', 'KDDC', 'KDFX',
        'KDGX', 'KDIX', 'KDLH', 'KDMX', 'KDTX', 'KDVN',
        'KEAX', 'KEMX', 'KENX', 'KEOX', 'KEPZ', 'KESX',
        'KEVX', 'KEWX', 'KEYX', 'KFCX', 'KFDR', 'KFDX',
        'KFFC', 'KFSD', 'KFSX', 'KFTG', 'KFWS', 'KGGW',
        'KGJX', 'KGLD', 'KGRB', 'KGRR', 'KGSP', 'KGWX',
        'KGYX', 'KHDX', 'KHGX', 'KHNX', 'KHPX', 'KHTX',
        'KICT', 'KICX', 'KILN', 'KILX', 'KIND', 'KINX',
        'KIWA', 'KIWX', 'KJAX', 'KJGX', 'KJKL', 'KLBB',
        'KLCH', 'KLGX', 'KLIX', 'KLNX', 'KLOT', 'KLRX',
        'KLSX', 'KLTX', 'KLVX', 'KLWX', 'KLZK', 'KMAF',
        'KMAX', 'KMBX', 'KMHX', 'KMKX', 'KMLB', 'KMOB',
        'KMPX', 'KMQT', 'KMRX']

start_date = datetime(2018, 9, 1)
end_date = datetime(2019, 9, 30)
s3 = boto3.client('s3', config=Config(signature_version=UNSIGNED))


sites_sliced = SITES[user_slice_start:user_slice_end] if user_slice_end else SITES

for site in sites_sliced:
    print(site)
    output_file = open("output/"+site+".asc", "a")
    today = start_date
    first = True
    while today <= end_date:
        date_str = today.strftime("%Y/%m/%d")
        print(date_str)
        prefix = f'{date_str}/{site}'
        response = s3.list_objects_v2(Bucket=NEXRAD_L2_BUCKET, Prefix=prefix)
        contents = response.get("Contents")
        if not contents:
            today = today + timedelta(days=1)
            continue
        first_today = True
        for item in contents:
            key = item.get("Key")
            print(key)
            if "MDM" in key:
                continue
            key_path = '/'.join(key.split('/')[:-1])
            if not os.path.exists(ROOT_PATH + key_path):
                os.makedirs(ROOT_PATH + key_path)
            if first:
                 s3.download_file(NEXRAD_L2_BUCKET, key, ROOT_PATH + key)
                 out = subprocess.run(["scansun_nexrad", ROOT_PATH + key], stdout=subprocess.PIPE)
                 lat = [float(item.split(':')[-1]) for item in str(out.stdout, 'utf-8').split('\n') if "Latitude" in item][0]
                 lon = [float(item.split(':')[-1]) for item in str(out.stdout, 'utf-8').split('\n') if "Longitude" in item][0]
                 first = False
            file_name = key.split('/')[-1]
            file_dt = datetime.strptime(file_name.replace('_', '').split('V06')[0][4:], "%Y%m%d%H%M%S")
            file_dt = file_dt.replace(tzinfo=tz.tzutc())
            if first_today:
                sun = Sun(lat, lon)
                today_sr = sun.get_sunrise_time(today)
                today_ss = sun.get_sunset_time(today)
                first_today = False
            sr_diff = (today_sr - file_dt).total_seconds()
            ss_diff = (today_ss - file_dt).total_seconds()
            if sr_diff < 10800 or ss_diff < 10800:
                if not os.path.exists(ROOT_PATH + key_path):
                    os.makedirs(ROOT_PATH + key_path)
                if not os.path.exists(ROOT_PATH + key):
                    s3.download_file(NEXRAD_L2_BUCKET, key, ROOT_PATH + key)
                out = subprocess.run(["scansun_nexrad", ROOT_PATH + key], stdout=subprocess.PIPE)
                out_str = str(out.stdout, 'utf-8')
                if "Sun" in out_str:
                    output_file.write(out_str)
                os.remove(ROOT_PATH + key)
        today = today + timedelta(days=1)
