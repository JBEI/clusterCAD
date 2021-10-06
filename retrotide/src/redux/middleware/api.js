// async actions and calls to the client / backend go here
// we need access to the client if any
// redux-thunk is imported toplevel so we should be able to access it and dispatch from here but we'll see

import axios from 'axios'; // this is the frontend client

const clusterCADSeqSearch = (molecule, token) => {
  console.log('hit api function search with ' + molecule);
  // return async function beginMoleculeSearch(dispatch, getState) {
    // const requestObject = { // object structure needed for clusterCAD request, probably shouldn't be here
      // csrfmiddlewaretoken: "",
      // sdf: molecule,
      // draw: 1,
      // cutoff: 0.0,
      // maxCompounds: 10,
    // };
    // make this a try/catch
    // const response = await client.post('/api', {integer: 0}); // callback
    // dispatch({ type: 'jobAdded', payload: {response} });
    axios.get('/api', {integer: 0})
    .then((response) => console.log(response))
    .catch((error) => console.log(error.config));
}

// export function synchWrapper (someInput) {
//   return async function beginMoleculeSearch(dispatch, getState) {
//     // example
//     const response = await client.post('/fakeApi/todos', { todo: initialTodo });
//     dispatch({ type: 'todos/todoAdded', payload: response.todo });
//   }
// }

export { clusterCADSeqSearch };
// then import this function in the component, which has no idea it's async
// remember to update reducers if dispatching an action