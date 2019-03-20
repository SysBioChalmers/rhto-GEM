#!/usr/bin/env bash

# Do NOT set -v or -x or your GitHub API token will be leaked!
set -ue # exit with nonzero exit code if anything fails

echo "Parse memote.ini for values."
deployment=$(awk -F '=' '{if (! ($0 ~ /^;/) && $0 ~ /deployment/) print $2}' memote.ini | tr -d ' ')
location=$(awk -F '=' '{if (! ($0 ~ /^;/) && $0 ~ /location/) print $2}' memote.ini | tr -d ' ')

echo "Configure Travis git user."
git config --global user.email "deploy@travis-ci.org"
git config --global user.name "Travis CI Deployment Bot"

# Generate the history report on the deployment branch.
echo "Generating new snapshot report."
memote report snapshot --filename "/tmp/index.html"
git checkout -- ./ComplementaryScripts/travis_deploy*.sh
git checkout "${deployment}"
mv "/tmp/index.html" ./

# Add, commit and push the files.
git add "index.html"
git commit -m "Travis report #${TRAVIS_BUILD_NUMBER}"
git push --quiet "https://${GITHUB_TOKEN}@github.com/${TRAVIS_REPO_SLUG}.git" "${deployment}" > /dev/null

echo "Your new report will be available at http://sysbiochalmers.github.io/rhto-GEM/ in a moment."
