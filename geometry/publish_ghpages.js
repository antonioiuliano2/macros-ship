/**
 * Given the following directory structure:
 *
 *   dist/
 *     index.html
 *     js/
 *       site.js
 *
 * The usage below will create a `gh-pages` branch that looks like this:
 *
 *   index.html
 *   js/
 *     site.js
 *
 */
var ghpages = require('gh-pages');
ghpages.publish('web',{branch: 'main', repo: 'git@github.com:antonioiuliano2/antonioiuliano2.github.io.git'}, function(err){});