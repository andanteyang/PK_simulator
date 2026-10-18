"""Optional real-browser regression checks; run against an existing static server."""
from playwright.sync_api import sync_playwright
from pathlib import Path
import csv
import math
import os
import tempfile

with sync_playwright() as pw, tempfile.TemporaryDirectory(prefix='pk-nca-check-') as output:
    browser = pw.chromium.launch(executable_path=os.environ.get('PK_CHROMIUM', '/usr/bin/chromium'), headless=True, args=['--no-sandbox'])
    page = browser.new_page(viewport={'width': 1440, 'height': 1100})
    errors = []
    page.on('pageerror', lambda error: errors.append(str(error)))
    page.goto(os.environ.get('PK_TEST_URL', 'http://127.0.0.1:8000/'), wait_until='networkidle')
    def settle():
        page.wait_for_timeout(180)
    def original_csv():
        with page.expect_download() as download:
            page.click('#export-csv')
        return Path(download.value.path()).read_text()
    for model in ['one', 'two']:
        for route in ['iv', 'po']:
            for dosing in ['single', 'multiple']:
                page.select_option('#model', model)
                page.select_option('#route', route)
                page.select_option('#dosing', dosing)
                settle()
                assert page.locator('#error').is_hidden()
                assert page.locator('#nca-panel').is_hidden()
                before = page.locator('#results').inner_text()
                chart = page.locator('#chart').evaluate('(c) => c.toDataURL()')
                curve_csv = original_csv()
                if dosing == 'multiple':
                    assert page.locator('#nca-enabled').is_disabled()
                    print('PASS NCA unavailable; original multiple-dose UI:', model, route)
                    continue
                page.check('#nca-enabled')
                assert page.locator('#nca-panel').is_visible()
                assert page.locator('#nca-error').is_hidden()
                assert page.locator('#results').inner_text() == before
                assert original_csv() == curve_csv
                assert page.locator('#nca-samples tbody tr').count() == 10
                page.click('#log')
                assert page.locator('#nca-error').is_hidden()
                page.click('#linear')
                page.uncheck('#nca-enabled')
                assert page.locator('#nca-panel').is_hidden()
                assert page.locator('#results').inner_text() == before
                assert page.locator('#chart').evaluate('(c) => c.toDataURL()') == chart
                print('PASS NCA toggle preserves results, original CSV and chart:', model, route)
    page.click('#reset')
    page.check('#nca-enabled')
    page.fill('#nca-times', '0, 2, 4, 8, 12, 24, 48')
    settle()
    assert 'No quantifiable IV 0 h postdose' in page.locator('#nca-content').inner_text()
    page.select_option('#nca-zero', 'post')
    assert '100,000' in page.locator('#nca-content').inner_text()
    assert page.locator('#nca-error').is_hidden()
    page.select_option('#nca-method', 'manual')
    assert 'At least three' in page.locator('#nca-content').inner_text()
    for t in [8, 12, 24]:
        page.check(f'[data-nca-time="{t}"]')
    assert 'Number of points: 3' in page.locator('#nca-content').inner_text()
    page.fill('#nca-lloq', '9000')
    settle()
    assert 'Manual points must' in page.locator('#nca-content').inner_text()
    assert 'BLQ' in page.locator('#nca-samples').inner_text()
    page.fill('#nca-lloq', '10')
    settle()
    page.select_option('#nca-method', 'auto')
    with page.expect_download() as download:
        page.click('#nca-export')
    records = list(csv.DictReader(Path(download.value.path()).read_text().splitlines()))
    assert len(records) == 7
    assert math.isclose(float(records[0]['auc0_inf_ng_h_ml']), 100000, rel_tol=1e-10)
    assert records[0]['Event state'] == 'post'
    assert 'lloq_ng_ml' in records[0]
    assert any(r['Terminal regression selection'] == 'selected' for r in records)
    page.fill('#nca-times', '0, 0')
    settle()
    assert 'must not repeat' in page.locator('#nca-error').inner_text()
    assert page.locator('#nca-export').is_disabled()
    assert page.locator('#error').is_hidden()
    page.fill('#nca-times', '0, 100')
    settle()
    assert 'outside simulation' in page.locator('#nca-error').inner_text()
    page.fill('#nca-times', '0, 2, 4, 8, 12, 24, 48')
    page.fill('#nca-lloq', '0')
    settle()
    assert 'LLOQ must be positive' in page.locator('#nca-error').inner_text()
    page.fill('#nca-lloq', '10')
    settle()
    # Quantified sample at 12 h, then the BLQ zero-time sample on a separate top rail.
    width = page.locator('#chart').evaluate('(c) => c.getBoundingClientRect().width')
    height = page.locator('#chart').evaluate('(c) => c.getBoundingClientRect().height')
    x12 = 69 + 12 / 72 * (width - 22 - 69)
    c12 = 10000 * math.exp(-1.2)
    y12 = height - 48 - c12 / 11300 * (height - 48 - 27)
    page.locator('#chart').hover(position={'x': x12, 'y': y12})
    assert 'Quantifiable' in page.locator('#tooltip').inner_text()
    page.select_option('#nca-zero', 'pre')
    page.locator('#chart').hover(position={'x': 69, 'y': 14})
    assert 'BLQ (< LLOQ; no reported concentration)' in page.locator('#tooltip').inner_text()
    page.click('#log')
    page.locator('#chart').hover(position={'x': 69, 'y': 14})
    assert 'BLQ' in page.locator('#tooltip').inner_text()
    with page.expect_download() as download:
        page.click('#export-png')
    assert Path(download.value.path()).read_bytes().startswith(b'\x89PNG')
    page.select_option('#dosing', 'multiple')
    assert page.locator('#nca-panel').is_hidden()
    page.select_option('#dosing', 'single')
    assert page.locator('#nca-panel').is_visible()
    page.select_option('#model', 'two')
    page.select_option('#route', 'po')
    page.set_viewport_size({'width': 390, 'height': 844})
    settle()
    assert page.evaluate('document.documentElement.scrollWidth <= innerWidth')
    page.screenshot(path='/tmp/pk-nca-mobile.png', full_page=True)
    page.set_viewport_size({'width': 1440, 'height': 1100})
    settle()
    page.screenshot(path='/tmp/pk-nca-desktop.png', full_page=True)
    assert not errors, errors
    print('PASS manual selection, IV zero state, diagnostics, NCA CSV, input errors, sample/BLQ tooltips, log scale, PNG, dosing transitions and mobile layout')
    browser.close()
