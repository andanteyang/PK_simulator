from playwright.sync_api import sync_playwright
from pathlib import Path
import os
with sync_playwright() as pw:
 browser=pw.chromium.launch(executable_path=os.environ.get('PK_CHROMIUM','/usr/bin/chromium'),headless=True,args=['--no-sandbox'])
 page=browser.new_page(viewport={'width':1440,'height':1100},device_scale_factor=1)
 errors=[]
 page.on('pageerror',lambda error:errors.append(str(error)))
 page.goto(os.environ.get('PK_TEST_URL','http://127.0.0.1:8000/'),wait_until='networkidle')
 for model in ['one','two']:
  for route in ['iv','po']:
   for dosing in ['single','multiple']:
    page.select_option('#model',model);page.select_option('#route',route);page.select_option('#dosing',dosing)
    page.wait_for_timeout(150)
    assert page.locator('#error').is_hidden(),page.locator('#error').inner_text()
    assert page.locator('.metric-value').count()>=3
    assert page.locator('#vc').count()==(model=='two')
    assert page.locator('#ka').count()==(route=='po')
    assert page.locator('#tau').count()==(dosing=='multiple')
    page.click('#log');assert page.locator('#log').get_attribute('aria-pressed')=='true'
    page.click('#linear')
    print('PASS UI',model,route,dosing)
 page.locator('#chart').hover(position={'x':300,'y':150})
 assert page.locator('#tooltip').is_visible()
 with page.expect_download() as info:page.click('#export-csv')
 download=info.value;download.save_as('/tmp/pk-simulation.csv')
 lines=Path('/tmp/pk-simulation.csv').read_text().splitlines()
 assert 'Event state' in lines[0] and len(lines)>600
 page.select_option('#route','iv')
 with page.expect_download() as info:page.click('#export-png')
 info.value.save_as('/tmp/pk-profile.png')
 assert Path('/tmp/pk-profile.png').read_bytes().startswith(b'\x89PNG')
 page.fill('#cl','0');page.wait_for_timeout(200)
 assert page.locator('#error').is_visible()
 assert page.locator('#export-csv').is_disabled()
 page.click('#reset');assert page.locator('#error').is_hidden()
 page.screenshot(path='/tmp/pk-desktop.png',full_page=True)
 page.set_viewport_size({'width':390,'height':844})
 page.wait_for_timeout(150)
 assert page.evaluate('document.documentElement.scrollWidth <= innerWidth')
 page.screenshot(path='/tmp/pk-mobile.png',full_page=True)
 assert not errors,errors
 print('PASS tooltip, CSV, PNG, invalid input, reset, mobile overflow, no browser errors')
 browser.close()
