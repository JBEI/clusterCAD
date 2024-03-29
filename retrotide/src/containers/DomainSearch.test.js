import React from 'react';
import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { Provider } from 'react-redux';
import { createStore, applyMiddleware, compose } from 'redux';
import thunk from 'redux-thunk';
import rootReducer from '../redux/reducers/rootReducer';
import DomainSearch from './DomainSearch';

// note that this is testing both DomainSearch and ModuleBuilder
// since they hand state back and forth

// DomainSearch is a connected component so we need to wrap
// in a Provider with a store. Redux docs suggest using the actual 
// store from the app

// should probably figure out how to render only once

const composeEnhancers = window.__REDUX_DEVTOOLS_EXTENSION_COMPOSE__ || compose;
const store = createStore(rootReducer, composeEnhancers(applyMiddleware(thunk)));

it('renders a loading and terminating module', () => {
  render(
    <Provider store={store}>
      <DomainSearch />
    </Provider>
  );

  const modules = screen.getAllByText(/Module [\d]/i);
  const loading = screen.getByText("loading");
  const terminating = screen.getByText("terminating");
  expect(modules).toHaveLength(2);
  expect(loading).toBeInTheDocument;
  expect(terminating).toBeInTheDocument;
});

it('adds new extending modules and deletes them', () => {
  render(
    <Provider store={store}>
      <DomainSearch />
    </Provider>
  );

  const addButton = screen.getByText(/Add Module +/);

  userEvent.click(addButton);
  userEvent.click(addButton);
  userEvent.click(addButton);

  var modules = screen.getAllByText(/Module [\d]/i);
  var extending = screen.getAllByText("extending");
  expect(modules).toHaveLength(5);
  expect(extending).toHaveLength(3);

  const deleteButton = screen.getAllByText(/X/)[0];

  userEvent.click(deleteButton);

  modules = screen.getAllByText(/Module [\d]/i);
  extending = screen.getAllByText("extending");
  expect(modules).toHaveLength(4);
  expect(extending).toHaveLength(2);
});

// submit button - ??

// tests for domain level

// expect to see 3 domain buttons in extending
// click button and expect other button to be disabled
// expect to see DH and KR in sandbox

// click KR in sandbox
// expect to see options buttons
// click options button
// expect to see clicked button selected
// click KR again
// expect to see modal close